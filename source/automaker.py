#! /usr/bin/python

class CppInterpreter:
    """
    Interprets c++ code
    """

    def __init__(self, cname, cflags):
        """
        Class constructor
        Input:
        cname - Compiler name, i.e. string used to invoke compiler
        cflag - Compilation flags
        """
        self.__cname = cname
        self.__cflags = cflags

    def SimpleCompile(self, sfname):
        """
        Creates an object file from a source file
        Input:
        sfname - Name of the source file
        """

        import os
        ef = os.system(self.__cname+' '+self.__cflags+' -c '+sfname)
        if ef!=0:
            raise NameError('Failed to compile')

    def NeedCompile(self,hfname):
        """
        Checkes whether a certain file needs to be compiled
        Returns true if so, and false if the file has already
        been compile (if a newer object file exists)
        """        

        res = 1
        if hfname[-4:]!='.hpp':
            print 'IfNeedCompile only accepts files with the extension .hpp'
            print 'Input file is '+hfname
            raise NameError('NotHeader')
        import os
        if os.path.exists(hfname.replace('.hpp','.o')):
            t_hpp = os.path.getmtime(hfname)
            t_o = os.path.getmtime(hfname.replace('.hpp','.o'))
            t_cpp = os.path.getmtime(hfname.replace('.hpp','.cpp'))
            if (t_o>t_cpp) and (t_o>t_hpp):
                res = 0
        return res                                   

    def RecursiveCompile(self, sfname):
        """
        Checks a source file for dependencies and compiles them recursively
        """

        f = open(sfname)
        for i in f:
            if i.find('#include')>-1 and i.find('"')>i.find('#define'):
                hname = i.split('"')[1]
                fname = hname.replace('hpp','cpp')
                if self.NeedCompile(hname) and sfname!=fname:
                    self.RecursiveCompile(fname)
                    self.RecursiveCompile(hname)
        f.close()
        if(sfname[-4:]=='.cpp'):
            self.SimpleCompile(sfname)

    def Link(self, sfname):
        """
        Links object files to an executable
        Input:
        sfname - Name of main file
        To - Do: The function assumes the input file has the extension .cpp
        At the moment I don't think that in this project I would need other
        extensions, so I decided to raise an error in case a file with a
        different extension is entered
        """

        if sfname.find('.cpp')==-1:
            raise NameError('automaker currently supports only .cpp files')

        import os
        ef = os.system(self.__cname+' '+self.__cflags+' -o '+
                       sfname.replace('.cpp','.exe')+' *.o')
        #ef = os.system(self.__cname+' -o '+sfname.replace('.cpp','.exe')+' *.o')
        if ef!=0:
            raise NameError('Failed to link '+sfname.replace('.cpp','.exe'))

    def Run(self, sfname):
        """
        Runs the executable
        This function assumes the executable was created using the Link function
        Input:
        sfname - Name of the main source file
        """

        if sfname.find('.cpp')==-1:
            raise NameError('automaker currently supports only .cpp files')

        import os
        os.system('./'+sfname.replace('.cpp','.exe'))

    def Interpret(self, sfname):
        """
        Compiles, links and executes
        Input:
        sfname - Name of the main source file
        """

        if sfname.find('.cpp')==-1:
            raise NameError('automaker currently supports only .cpp files')

        self.RecursiveCompile(sfname)
        self.Link(sfname)
        self.Run(sfname)

import sys
if __name__=="__main__":
    fname = sys.argv[1]

    import os
    os.system('rm *.o *.exe')
    myinterp = CppInterpreter('g++','')
    myinterp.Interpret(fname)

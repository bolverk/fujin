#! /usr/bin/python

class Tester:
    def __init__(self, tfname, srcdir):
        """
        Class constructor
        Input:
        tfname - Name of temporary folder
        srcdir - Path to source
        """

        self.__tfname = tfname
        self.__srcdir = srcdir
        return

    def RunSingleTest(self, dname):
        """
        Runs a single test
        Input:
        dname - Path to directory
        """

        import os
        import imp
        import re
        temp_dir = self.__tfname+'_'+dname.replace('.','').replace('/','_')
        temp_dir = re.sub(r'_$','',temp_dir)
        os.system('touch '+temp_dir)
        os.system('rm -r '+temp_dir)
        os.mkdir(temp_dir)
        os.system('cp '+self.__srcdir+'/* '+temp_dir)
        os.system('cp '+dname+'/* '+temp_dir)
        home_dir = os.getcwd()
        os.chdir(temp_dir)
        automaker = imp.load_source('automaker',os.getcwd()+'/automaker.py')
        intrp = automaker.CppInterpreter(\
            'g++', '-Wall -Wextra -pedantic -Werror -pedantic -ansi -Wfatal-errors -O2')
        intrp.Interpret('test.cpp')
        mod = imp.load_source('mod',os.getcwd()+'/test.py')
        res = mod.main()
        os.chdir(home_dir)
        return res
    
    def DelTempDirs(self):
        """
        Deletes the temporary folder and all files within
        """

        import os
        os.system('touch '+self.__tfname)
        os.system('rm -r '+self.__tfname+'*')
        return

def ListLeaves(dname):
    """
    Returns a list with the names of all the leaves
    (folders without subfolders)
    """

    import os
    res = []
    for roots, dirs, files in os.walk(dname):
        if roots.find('.svn')==-1 and len(roots)>1:
            if os.path.isfile(roots+'/test.py'):
                res.append(roots)
    return res

import sys
import os
import multiprocessing

tst = Tester('temp','../source')

def run_test_error_handle(tester,test_dir):
    home_dir = os.getcwd()
    try:
        tf = tst.RunSingleTest(test_dir)
        if tf:
            print test_dir+' ... passed'
        else:
            print test_dir+' ... ran and failed'
    except Exception, err:
        print test_dir+' ... failed to run'
        print err
        os.chdir(home_dir)

if len(sys.argv)==1:
    tst.DelTempDirs()
    ll = ListLeaves('.')
    my_pool = multiprocessing.Pool()
    for leaf in ll:
        my_pool.apply_async(run_test_error_handle,(tst,leaf))
        #run_test_error_handle(tst,leaf)
    my_pool.close()
    my_pool.join()
    tst.DelTempDirs()
elif len(sys.argv)==2:
    run_test_error_handle(tst,sys.argv[1])
else:
    raise NameError('Unknown option')

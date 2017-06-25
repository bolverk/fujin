/*! \mainpage Fujin - Relativistic Hydrodynamic Simulation
 *
 * \htmlonly
 * <P>
 * <img src="../../source/fujin.jpg" alt="Fujin statue" align="left" hspace="10"/>
 * </P>
 * \endhtmlonly
 *
 * \section intro_sec Introduction
 *
 * Fujin is a Lagrangian, one dimensional, Godunov method special relativistic hydrodynamic simulation, written in c++.
 *
 * \section code_standard Coding Standards
 *
 * \subsection style Coding Style
 *
 * I've arbitrarily chosen the coding style described <a href url="http://geosoft.no/development/cppstyle.html">here</a>.
 *
 * \subsection warnings Compilation Warnings
 *
 * All C++ code in the Fujin project should compile with <a href url="http://gcc.gnu.org/">g++</a>, using the following flags\n
 * \code 
 * -Wall -Wextra -pedantic -Werror -ansi 
 * \endcode 
 *
 * \subsection documentation Documentation 
 *
 * Everything must be documented with <a href url="http://www.stack.nl/~dimitri/doxygen/"> Doxygen </a>. Every function must contain a description of the input and the output. Any additional information (bugs, to-dos etc) is encouraged.
 * In fact, when it comes to documentation - the more the merrier.\n
 * Mathematica files should be documented according to the instruction <a href url="http://reference.wolfram.com/mathematica/tutorial/DocumentationConstructs.html"> here </a>.\n
 * Python files should be documented using <a href url="http://www.python.org/dev/peps/pep-0257/"> docstring </a>.
 *
 * \subsection autotest Tests
 * For every major function, the developer should write a program that tests whether the function works properly. Whether a function is complex enough to warrant a test, and the nature of the test is at the developer's discretion.
 * Whenever the developer encouters and fixes a bug, she should write a test that checks for that bug. This "bug trap" would prevent that same bug from recurring.
 * As for the technical details of adding a test, a test is composed of two files, test.cpp and test.py. The first file does some calculation with the function we wish to test, and writes the result to the file "res.txt". The file "test.py" verifies that the result is correct.
 * These two files should be in a folder with the name of the tested function. That subfolder should be in a folder with the name of the cpp file. All the tests should be in the "tests" folder.
 * To run all test, cd to the test folder and type \n
 * \code
 * ./run_all_tests.py
 * \endcode
 * This will create a new folder "temp", copy the files there, compile, link, run the program and compare the results. If a test went well, it will echo the path to the test and the string 'passed'. If something went wrong, it will also echo the path, but add the string 'fail'. 
 * In order to debug a certain test, say, the one in the folder igrs/hugoniotvelocity, use the command \n
 * \code
 * ./run_all_test.py debug igrs/hugoniotvelocty
 * \endcode
 * This will leave everything in the folder "temp". When you debug, remember that the files inside "temp" are copied from other folders, so changes made to them will be erased (correcting the files there won't solve the problem). Delete the "temp" folder when you're done.
 *
 * \subsection versioncontrol Version Control
 *
 * Once a patch is written according to the style, checked for compiler errors, documented and tested, it should be commited to the repository using the <a href url="http://subversion.tigris.org/"> Subversion </a> version control system.
 *
 * \section misc Miscellaneous
 *
 * \subsection background Background
 *
 * My name is Almog Yalinewich, and I've prepared this code as a part of my doctoral research, under the supervision of Prof. Re'em Sari. This code is largely based on the code described in <a href url="adsabs.harvard.edu/abs/1999ApJ...513..669k"> S. Kobayashi, T. Piran, R. Sari, "Hydrodynamics of a Relativistic Fireball: The Complete Evolution", Astrophys. Jour. 513, 669-678 (1999) </a>.
 *
 * \subsection namesake Namesake
 *
 * <a href url="http://en.wikipedia.org/wiki/F%C5%ABjin">Fujin</a> is the Japanese god of wind and storms. At some point I thought this simulation will be used to describe shock waves in stellar wind.
 *
 */

import glob
import os
from argparse import ArgumentParser

debug = ARGUMENTS.get('debug',0)
source = ARGUMENTS.get('source',None)
target = ARGUMENTS.get('target',None)
compiler = ARGUMENTS.get('compiler','g++')

if compiler=='g++':
    cflags = '-Wfatal-errors'
    if int(debug):
        cflags +=' -O0 -g -pg -DDEBUG_MODE'
    else:
        cflags+=' -O3'
elif compiler=='clang++':
    cflags = '-Weverything -Werror -ferror-limit=1 -Wno-error=padded'
    if int(debug):
        cflags += ' -O0 -g -pg -DDEBUG_MODE'
    else:
        cflags += ' -O3 -march=native'
env = Environment(ENV = os.environ,
                  CXX=compiler,
                  CPPPATH=os.environ['FUJIN_ROOT']+'/source',
                  LIBPATH=['.',os.environ['HDF5_LIB_PATH']],
                  LIBS=['fujin','hdf5','hdf5_cpp'],
                  CXXFLAGS=cflags)
lib_file = env.Library('fujin',glob.glob(os.environ['FUJIN_ROOT']+'/source/*.cpp'))
if None!=source and None!=target:
    env.Program(target,[source])

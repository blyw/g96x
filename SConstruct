import os
import socket
import getpass
import time

build_dir = time.strftime("%Y%m%d%H%M%S", time.gmtime()) + '_' + getpass.getuser() + '_' + socket.gethostname()

env = Environment(ENV = {'PATH' : os.environ['PATH'], 'LD_LIBRARY_PATH' : os.environ['LD_LIBRARY_PATH']})

env.Append(CCFLAGS=['-O3', '-std=c++11', '-march=corei7-avx'])
env.Append(CPPPATH=['.',
 '/pool/lai/prog/shared-lib/eigen-3.2.0',
 '/pool/lai/prog/ubuntu-12.04/libxml2-2.9.1/include/libxml2',
 ])
env.Append(LIBPATH=['/pool/lai/prog/ubuntu-12.04/libxml2-2.9.1/lib','/pool/lai/prog/ubuntu-12.04/zlib-1.2.8/lib'])

env.Append(LINKFLASGS=['-lz', '-lxml2', '-pthread'])
env.Append(LIBS = ['libxml2','pthread','z'])

Export('env')
SConscript('frameoutX/SConscript', variant_dir='.build_' + build_dir)


#!/bin/bash
echo have to set PI_DIR and PYTHON_LIB_DIR and SWIG_BASE_DIR
if [ -n $SWIG_BASE_DIR ]; then
    echo SWIG_BASE_DIR default is /usr/local
    export SWIG_BASE_DIR=/usr/local
else
    echo SWIG_BASE_DIR is $SWIG_BASE_DIR
fi
if [ -n $PYTHON_INCLUDE_DIR ]; then
    export PYTHON_INCLUDE_DIR=/usr/include/python2.3
    echo PYTHON_INCLUDE_DIR default is $PYTHON_INCLUDE_DIR
else
    echo PYTHON_INCLUDE_DIR is $PYTHON_INCLUDE_DIR
fi
if [ -n $PYTHON_LIB_DIR ]; then
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PYTHON_LIB_DIR}
else
    echo Please set PYTHON_LIB_DIR
fi
export PYTHONVERSION=2.3
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${G4WORKDIR}/tmp/${G4SYSTEM}/TiaraWrapper:${G4WORKDIR}/tmp/${G4SYSTEM}/G4KernelWrapper:${G4WORKDIR}/tmp/${G4SYSTEM}/CLHEPWrapper
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${G4WORKDIR}/tmp/${G4SYSTEM}/tiara
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${G4WORKDIR}/lib/${G4SYSTEM}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PYTHON_LIB_DIR}
export TIARA_BASE=`pwd`
export TIARASCRIPTS=${TIARA_BASE}/source/py_modules
export PYTHONPATH=${PYTHONPATH}:${LD_LIBRARY_PATH}
export PYTHONPATH=${PYTHONPATH}:${TIARA_BASE}/run:${TIARASCRIPTS}
export PYTHONPATH=${PYTHONPATH}:${TIARA_BASE}/source/TiaraWrapper
export PYTHONPATH=${PYTHONPATH}:${TIARA_BASE}/source/G4KernelWrapper
export PYTHONPATH=${PYTHONPATH}:${TIARA_BASE}/source/CLHEPWrapper
if [ -n $SWIG_VERSION ]; then
    export SWIG_VERSION=1.3.15
    echo setting SWIG_VERSION to default $SWIG_VERSION
else
    echo SWIG_VERSION is $SWIG_VERSION
fi
export SWIG_INCDIRS="-I${SWIG_BASE_DIR}/lib/swig_lib -I${SWIG_BASE_DIR}/lib/swig_lib/python"
if [ -n $SWIG ]; then
    export SWIG=${SWIG_BASE_DIR}/bin/swig
    echo setting SWIG executable to $SWIG
else
    echo SWIG executable is $SWIG
fi
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${G4WORKDIR}/tmp/${G4SYSTEM}/CLHEPWrapper
export PYTHONPATH=${PYTHONPATH}:${TIARA_BASE}/source/CLHEPWrapper
echo if python2.3 is not picked up then execute following command
echo "'cp $PI_DIR/lib/libpython2.3.so $G4LIB/$G4SYSTEM/.'"
echo finished envCommon.sh

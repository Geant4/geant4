#!/bin/bash
echo have to set PI_DIR and PYTHON_LIB_DIR
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
export SWIG_VERSION=1.3.15
export SWIG_INCDIRS="-I${SWIG_BASE_DIR}/lib/swig_lib -I${SWIG_BASE_DIR}/lib/swig_lib/python"
export SWIG=${SWIG_BASE_DIR}/bin/swig
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${G4WORKDIR}/tmp/${G4SYSTEM}/CLHEPWrapper
export PYTHONPATH=${PYTHONPATH}:${TIARA_BASE}/source/CLHEPWrapper
echo if python2.3 is not picked up then execute following command
echo "'cp $PI_DIR/libpython2.3.so $G4LIB/$G4SYSTEM/.'"
echo finished envCommon.sh
#!/bin/tcsh -x

#set verbose = 0

if (${?PI_BASE_DIR} == 1) then
   setenv PI_VERSION 1_2_1
   setenv PI_VER 1.2.1
   setenv PLATF rh73_gcc32
   setenv PATH ${PATH}:${PI_BASE_DIR}/${PI_VER}/app/releases/PI/PI_${PI_VER}/${PLATF}/bin
   setenv PATH ${PI_BASE_DIR}/Linux-g++/bin:$PATH
   eval `aida-config -r csh`
#
# python from PI
#
  setenv PYTHONVERSION 2.2
  setenv PYTHON_INCLUDE_DIR ${PI_BASE_DIR}/${PI_VER}/external/Python/2.2.2/${PLATF}/include/python${PYTHONVERSION}
  setenv PYTHON_LIB_DIR ${PI_BASE_DIR}/${PI_VER}/external/Python/2.2.2/${PLATF}/lib/python${PYTHONVERSION}/config
#
# python for PI
if (${?PYTHONPATH} == 0) then
  setenv PYTHONPATH
endif
  setenv PYTHONPATH ${PI_BASE_DIR}/Linux-g++/python:$PYTHONPATH
  setenv PYTHONPATH ${PI_BASE_DIR}/${PI_VER}/app/releases/SEAL/SEAL_1_3_4/${PLATF}/lib:$PYTHONPATH
#
#
else   # not using PI

  echo "-- WARNING: histograms are not activated !"
  echo "            Either PI_BASE_DIR are not set."
  #
  if (   ${?PYTHONVERSION} == 1 && \
         ${?PYTHON_BASE_DIR} == 1 && \
         ${?CLHEP_BASE_DIR} == 1) then  # settings without PI
    setenv PYTHON_LIB_DIR ${PYTHON_BASE_DIR}/lib/python${PYTHONVERSION}/config
    setenv PYTHON_INCLUDE_DIR ${PYTHON_BASE_DIR}/include/python${PYTHONVERSION}
    #
    if ( ! -d $PYTHON_LIB_DIR ) then
      echo -- ERROR: no pyhton lib/config directory: $PYTHON_LIB_DIR
    endif
    if ( ! -d $PYTHON_INCLUDE_DIR ) then
      echo -- ERROR: no pyhton include directory: $PYTHON_INCLUDE_DIR
    endif
  else    # environment not completed in case no PI is used
    echo "-- ERROR: PYTHONVERSION or PYTHON_BASE_DIR or CLHEP_BASE_DIR not set!"    exit
  endif

endif  # end settings for not using PI

if ( ${?LD_LIBRARY_PATH} == 0 ) then
  setenv LD_LIBRARY_PATH
endif

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${G4WORKDIR}/lib/${G4SYSTEM}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${PYTHON_LIB_DIR}
setenv LD_LIBRARY_PATH ${PI_BASE_DIR}/${PI_VER}/app/releases/SEAL/SEAL_1_3_4/rh73_gcc32/lib:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH /afs/cern.ch/sw/lcg/external/Python/2.2.2/rh73_gcc32/lib:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH ${PI_BASE_DIR}/${PI_VER}/external/Boost/1.30.2/rh73_gcc32/lib:$LD_LIBRARY_PATH

if (${?PYTHONPATH} == 0) then
  setenv PYTHONPATH
endif
setenv PYTHONPATH ${PYTHONPATH}:${LD_LIBRARY_PATH}

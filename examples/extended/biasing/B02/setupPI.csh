#!/bin/tcsh -x

#set verbose = 0

if (${?PI_BASE_DIR} == 1) then
   setenv PI_VER 1.2.6-lite
   setenv PLATF slc3_gcc323 
   setenv PI_DIR ${PI_BASE_DIR}/${PI_VER}/${PLATF}
   setenv PATH ${PI_DIR}/bin:${PATH}
   eval `aida-config -r csh`
#
# python from PI
#
  setenv PYTHON_LIB_DIR ${PI_DIR}/lib
#
# python for PI
if (${?PYTHONPATH} == 0) then
  setenv PYTHONPATH
endif
  setenv PYTHONPATH ${PI_DIR}/python:$PYTHONPATH
  setenv PYTHONPATH ${PI_DIR}/lib:$PYTHONPATH
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
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${G4WORKDIR}/tmp/${G4SYSTEM}}/exampleB02
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${PYTHON_LIB_DIR}
setenv LD_LIBRARY_PATH ${PI_DIR}/lib:$LD_LIBRARY_PATH

if (${?PYTHONPATH} == 0) then
  setenv PYTHONPATH
endif
#setenv PYTHONPATH ${PYTHONPATH}:${LD_LIBRARY_PATH}

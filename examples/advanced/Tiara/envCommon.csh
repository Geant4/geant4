# $Id: envCommon.csh,v 1.13 2005/03/17 19:48:27 daquinog Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-09-00 $
# -------------------------------------------------------------------
# Before sourcing this script make sure you have set the 
# environment variables according to the description in README.
# -------------------------------------------------------------------

# setup in case PI is used --------------------------------------

if (${?PI_BASE_DIR} == 1 && ${?SWIG_BASE_DIR} == 1) then

  setenv PI_VERSION 1_3_0
  setenv PI_VER 1.3.0
  setenv PI_ARCH slc3_ia32_gcc323
  setenv PYTHONVERSION 2.3
  setenv PI_DIR ${PI_BASE_DIR}/${PI_VER}/app/releases/PI/PI_${PI_VERSION}/${PI_ARCH}
  setenv PATH ${PATH}:${PI_DIR}/bin
  setenv PATH ${PI_BASE_DIR}/${PI_VER}/external/Python/${PYTHONVERSION}.4/${PI_ARCH}/bin:$PATH
  eval `aida-config -r csh`
  #
  # python from PI
  #
  setenv PYTHON_INCLUDE_DIR ${PI_BASE_DIR}/${PI_VER}/external/Python/${PYTHONVERSION}.4/${PI_ARCH}/include/python${PYTHONVERSION}
  setenv PYTHON_LIB_DIR ${PI_BASE_DIR}/${PI_VER}/external/Python/${PYTHONVERSION}.4/${PI_ARCH}/lib/python${PYTHONVERSION}/config
  #
  # python for PI
if (${?PYTHONPATH} == 0) then
  setenv PYTHONPATH
endif
  setenv PYTHONPATH ${PI_DIR}/python:$PYTHONPATH
  setenv PYTHONPATH ${PI_BASE_DIR}/${PI_VER}/app/releases/SEAL/SEAL_1_6_0/${PI_ARCH}/python:$PYTHONPATH
  #
  # setup the swig command
  #
  setenv SWIG_VERSION 1.3.15
  setenv SWIG_INCDIRS "-I${SWIG_BASE_DIR}/lib/swig-${SWIG_VERSION} -I${SWIG_BASE_DIR}/lib/swig-${SWIG_VERSION}/python"
  setenv SWIG ${SWIG_BASE_DIR}/bin/swig-${SWIG_VERSION}

else   # not using PI

  echo "-- WARNING: histograms are not activated !"
  echo "            Either PI_BASE_DIR or SWIG_BASE_DIR are not set."
  #
  if (   ${?PYTHONVERSION} == 1 && \
         ${?PYTHON_BASE_DIR} == 1 && \
         ${?SWIG_BASE_DIR} == 1 && \
         ${?SWIG_VERSION} == 1 && \
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
    #
    setenv SWIG_INCDIRS "-I${SWIG_BASE_DIR}/lib/swig-${SWIG_VERSION} -I${SWIG_BASE_DIR}/lib/swig-${SWIG_VERSION}/python"
    if ( -x  ${SWIG_BASE_DIR}/bin/swig-${SWIG_VERSION} ) then
      setenv SWIG ${SWIG_BASE_DIR}/bin/swig-${SWIG_VERSION}
    else if ( -x ${SWIG_BASE_DIR}/bin/swig ) then
      setenv SWIG ${SWIG_BASE_DIR}/bin/swig
    else 
      echo -- ERROR: could not find swig executable !
    endif
  else    # environment not completed in case no PI is used
    echo "-- ERROR: PYTHONVERSION or PYTHON_BASE_DIR or SWIG_BASE_DIR or SWIG_VERSION or CLHEP_BASE_DIR not set!"
    exit
  endif

endif    # end settings for not using PI

#
# common settings 
#

if (${?LD_LIBRARY_PATH} == 0) then
  setenv LD_LIBRARY_PATH
endif
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${G4WORKDIR}/tmp/${G4SYSTEM}/TiaraWrapper
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${G4WORKDIR}/tmp/${G4SYSTEM}/G4KernelWrapper
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${G4WORKDIR}/tmp/${G4SYSTEM}/CLHEPWrapper
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${G4WORKDIR}/tmp/${G4SYSTEM}/tiara
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${G4WORKDIR}/lib/${G4SYSTEM}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${PYTHON_LIB_DIR}
#setenv LD_LIBRARY_PATH ${PI_BASE_DIR}/${PI_VER}/app/releases/SEAL/SEAL_1_3_4/rh73_gcc32/lib:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH ${PI_BASE_DIR}/${PI_VER}/external/Python/${PYTHONVERSION}.4/${PI_ARCH}/lib:${LD_LIBRARY_PATH}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${PI_BASE_DIR}/${PI_VER}/external/Boost/1.31.0/${PI_ARCH}/lib
setenv TIARA_BASE `pwd`
setenv TIARASCRIPTS ${TIARA_BASE}/source/py_modules
if (${?PYTHONPATH} == 0) then
  setenv PYTHONPATH
endif
setenv PYTHONPATH ${PYTHONPATH}:${LD_LIBRARY_PATH}
setenv PYTHONPATH ${PYTHONPATH}:${TIARA_BASE}/run:${TIARASCRIPTS}
setenv PYTHONPATH ${PYTHONPATH}:${TIARA_BASE}/source/TiaraWrapper
setenv PYTHONPATH ${PYTHONPATH}:${TIARA_BASE}/source/G4KernelWrapper
setenv PYTHONPATH ${PYTHONPATH}:${TIARA_BASE}/source/CLHEPWrapper
#

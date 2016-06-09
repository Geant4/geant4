# $Id: envCommon.csh,v 1.9 2004/06/09 15:04:34 daquinog Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-06-02 $
# -------------------------------------------------------------------
# Before sourcing this script make sure you have set the 
# environment variables according to the description in README.
# -------------------------------------------------------------------

# setup in case Anaphe is used --------------------------------------

if (${?PI_BASE_DIR} == 1 && ${?SWIG_BASE_DIR} == 1) then

  setenv PI_VERSION 1_2_1
  setenv PI_VER 1.2.1
  setenv PI_ARCH rh73_gcc32
  setenv PATH ${PATH}:${PI_BASE_DIR}/${PI_VER}/app/releases/PI/PI_${PI_VER}/${PI_ARCH}/bin
  setenv PATH ${PI_BASE_DIR}/Linux-g++/bin:$PATH
  eval `aida-config -r csh`
  #
  # python from PI
  #
  setenv PYTHONVERSION 2.2
  setenv PYTHON_INCLUDE_DIR ${PI_BASE_DIR}/${PI_VER}/external/Python/2.2.2/${PI_ARCH}/include/python${PYTHONVERSION}
  setenv PYTHON_LIB_DIR ${PI_BASE_DIR}/${PI_VER}/external/Python/2.2.2/${PI_ARCH}/lib/python${PYTHONVERSION}/config
  #
  # python for PI
if (${?PYTHONPATH} == 0) then
  setenv PYTHONPATH
endif
  setenv PYTHONPATH ${PI_BASE_DIR}/Linux-g++/python:$PYTHONPATH
  setenv PYTHONPATH ${PI_BASE_DIR}/${PI_VER}/app/releases/SEAL/SEAL_1_3_4/rh73_gcc32/lib:$PYTHONPATH
  #
  # setup the swig command
  #
  setenv SWIG_VERSION 1.3.15
  setenv SWIG_INCDIRS "-I${SWIG_BASE_DIR}/lib/swig-${SWIG_VERSION} -I${SWIG_BASE_DIR}/lib/swig-${SWIG_VERSION}/python"
  setenv SWIG ${SWIG_BASE_DIR}/bin/swig-${SWIG_VERSION}

else   # not using Anaphe

  echo "-- WARNING: histograms are not activated !"
  echo "            Either PI_BASE_DIR or SWIG_BASE_DIR are not set."
  #
  if (   ${?PYTHONVERSION} == 1 && \
         ${?PYTHON_BASE_DIR} == 1 && \
         ${?SWIG_BASE_DIR} == 1 && \
         ${?SWIG_VERSION} == 1 && \
         ${?CLHEP_BASE_DIR} == 1) then  # settings without Anaphe
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
  else    # environment not completed in case no Anaphe is used
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
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${G4WORKDIR}/tmp/${G4SYSTEM}/tiaraPhysicsPackaging
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${G4WORKDIR}/tmp/${G4SYSTEM}/tiaraPhysicsLists
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${G4WORKDIR}/lib/${G4SYSTEM}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${PYTHON_LIB_DIR}
setenv LD_LIBRARY_PATH ${PI_BASE_DIR}/${PI_VER}/app/releases/SEAL/SEAL_1_3_4/rh73_gcc32/lib:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH /afs/cern.ch/sw/lcg/external/Python/2.2.2/rh73_gcc32/lib:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH ${PI_BASE_DIR}/${PI_VER}/external/Boost/1.30.2/rh73_gcc32/lib:$LD_LIBRARY_PATH
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

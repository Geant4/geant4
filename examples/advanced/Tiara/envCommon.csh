# $Id: envCommon.csh,v 1.8 2003/12/08 17:53:26 gcosmo Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-06-00 $
# -------------------------------------------------------------------
# Before sourcing this script make sure you have set the 
# environment variables according to the description in README.
# -------------------------------------------------------------------

# setup in case Anaphe is used --------------------------------------

if (${?ANAPHE_SCRIPTS} == 1 && ${?ANAPHETOP} == 1 && ${?AIDA_DIR} == 1 ) then

  if ( -f $ANAPHE_SCRIPTS/setupAnaphe.csh ) then 
    source $ANAPHE_SCRIPTS/setupAnaphe.csh 
  else if ( -f $ANAPHE_SCRIPTS/setup-Anaphe.csh) then
    source $ANAPHE_SCRIPTS/setup-Anaphe.csh
  else
    echo "envCommon.csh: ERROR: could not find a setup file for anaphe"
    exit
  endif
  setenv PATH ${ANAPHE_SCRIPTS}:$PATH
  setenv ANAPHESPECDIR `(cd $ANAPHE_REL_DIR/../; pwd)`
  #
  # python from Anaphe
  #
  setenv PYTHONVERSION 2.2
  setenv PYTHON_INCLUDE_DIR ${ANAPHESPECDIR}/PublicDomainPackages/2.0.0/include/python${PYTHONVERSION}
  setenv PYTHON_LIB_DIR ${ANAPHESPECDIR}/PublicDomainPackages/2.0.0/lib/python${PYTHONVERSION}/config
  #
  # setup the swig command
  #
  setenv SWIG_VERSION 1.3.15
  setenv SWIG_BASE_DIR ${ANAPHESPECDIR}/PublicDomainPackages/2.0.0
  setenv SWIG_INCDIRS "-I${SWIG_BASE_DIR}/lib/swig-${SWIG_VERSION} -I${SWIG_BASE_DIR}/lib/swig-${SWIG_VERSION}/python"
  setenv SWIG ${SWIG_BASE_DIR}/bin/swig-${SWIG_VERSION}
  #
  # CLHEP from Anaphe (shared libraries are required)
  #
  setenv CLHEP_BASE_DIR ${ANAPHETOP}/specific/${PLATF}/CLHEP/1.8.0.0
  #
  # add scripts path to PATH
  #
  setenv PATH ${ANAPHETOP}/share/LHCXX/${ANAPHE_VERSION}/scripts:${PATH}

else   # not using Anaphe

  echo "-- WARNING: histograms are not activated !"
  echo "            Either ANAPHE_SCRIPTS, ANAPHETOP or AIDA_DIR are not set."
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

endif    # end settings for not using Anaphe

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

# $Id: envCommon.csh,v 1.3 2003-06-16 15:49:06 dressel Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: not supported by cvs2svn $
# -------------------------------------------------------------------

# before sourcing this script make sure you have set the 
# environment variables according to the description in README


# in case Anaphe is used via afs.
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

# python from Anaphe
setenv PYTHONVERSION 2.2

setenv PYTHON_INCLUDE_DIR ${ANAPHESPECDIR}/PublicDomainPackages/2.0.0/include/python${PYTHONVERSION}
setenv PYTHON_LIB_DIR ${ANAPHESPECDIR}/PublicDomainPackages/2.0.0/lib/python${PYTHONVERSION}/config

# set the swig command
setenv SWIG ${ANAPHESPECDIR}/PublicDomainPackages/2.0.0/bin/swig-1.3.15

# CLHEP access
setenv CLHEP_BASE_DIR ${ANAPHETOP}/specific/${PLATF}/CLHEP/1.8.0.0/

else

echo "envCommon.csh: INFO: Anaphe is not used since ANAPHE_SCRIPTS or ANAPHETOP or AIDA_DIR is not set"

endif

if (   ${?PYTHONVERSION} == 1 && \
       ${?PYTHON_INCLUDE_DIR} == 1 && \
       ${?PYTHON_LIB_DIR} == 1 && \
       ${?SWIG} == 1 && \
       ${?CLHEP_BASE_DIR} == 1) then

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

if (${?ANAPHETOP} == 1) then 
setenv PATH ${ANAPHETOP}/share/LHCXX/${ANAPHE_VERSION}/scripts:${PATH}
endif

else 

echo "envCommon.csh: ERROR: PYTHONVERSION or PYTHON_INCLUDE_DIR or PYTHON_LIB_DIR or SWIG or CLHEP_BASE_DIR  not set!"

endif

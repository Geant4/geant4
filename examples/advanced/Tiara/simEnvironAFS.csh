# set the variables related to ANAPHE
setenv ANAPHE_SCRIPTS /afs/cern.ch/sw/lhcxx/share/LHCXX/5.0.4/scripts
source $ANAPHE_SCRIPTS/setupAnaphe.csh  #setup-Anaphe.csh in later vers.
setenv PATH ${ANAPHE_SCRIPTS}:$PATH
setenv ANAPHESPECDIR `(cd $ANAPHE_REL_DIR/../; pwd)`



# set variables related to python
setenv PYTHONVERSION 2.2

setenv PYTHON_INCLUDE_DIR ${ANAPHESPECDIR}/PublicDomainPackages/2.0.0/include/python${PYTHONVERSION}
setenv PYTHON_LIB_DIR ${ANAPHESPECDIR}/PublicDomainPackages/2.0.0/lib/python${PYTHONVERSION}/config


# this part should be fine in most cases
setenv G4LIB_BUILD_SHARED 1
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




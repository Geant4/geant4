#!/bin/tcsh -f

setenv APPLICATIONDIR /grid/app/geant4_validation
setenv G4WORKDIR /grid/app/geant4_validation/geant4
source ${APPLICATIONDIR}/env.csh 
if ( ${?PATH} ) then 
    setenv PATH ${PATH}:/bin:/usr/bin
    echo $PATH
else
    setenv PATH /bin:/usr/bin
    echo $PATH
endif
setenv ROOTSYS /grid/app/geant4_validation/root
#setenv PATH ${ROOTSYS}/bin:${PATH}
setenv LD_LIBRARY_PATH ${ROOTSYS}/lib:${LD_LIBRARY_PATH}
echo ${LD_LIBRARY_PATH} 
setenv G4EXE ${G4WORKDIR}/bin/${G4SYSTEM}
setenv ROOTSCRIPTS  ${G4WORKDIR}/tests/test47

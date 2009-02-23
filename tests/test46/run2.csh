#/bin/csh

cd $VFEM/test46
mkdir -p $REFERENCE
cd $REFERENCE

setenv PRIMARYBEAM proton

setenv PHYSLIST QBBC
#source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QBBCG
#source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QBBCF
#source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QGSP_BERT
#source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QGSP_BERT_EMV
source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST FTFP_BERT
source $G4INSTALL/tests/test46/run.csh

#

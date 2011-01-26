#/bin/csh

cd $VFHAD/test46
mkdir -p $REFERENCE
cd $REFERENCE

setenv PRIMARYBEAM proton
mkdir -p $PRIMARYBEAM
cd $PRIMARYBEAM

setenv PHYSLIST QGSP_BERT
source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QGSP_BERT_EMV
source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST FTFP_BERT
source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST FTFP_BERT_EMV
source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QGSP_FTFP_BERT
source $G4INSTALL/tests/test46/run.csh

cd ../
setenv PRIMARYBEAM pi-
mkdir -p $PRIMARYBEAM
cd $PRIMARYBEAM

setenv PHYSLIST QGSP_BERT
source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QGSP_BERT_EMV
source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST FTFP_BERT
source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST FTFP_BERT_EMV
source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QGSP_FTFP_BERT
source $G4INSTALL/tests/test46/run.csh

#

#/bin/csh

cd $VFEM/test46
mkdir -p $REFERENCE
cd $REFERENCE

setenv PRIMARYBEAM pi-
mkdir -p $PRIMARYBEAM
cd $PRIMARYBEAM
 
setenv PHYSLIST QBBC
#source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QBBCG
#source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QBBCF
#source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QGSP_BERT_EML
<<<<<<< run1.csh
<<<<<<< run1.csh
source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QGSP_BERT_EMLSN
#source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QBBC_XGGSN
#source $G4INSTALL/tests/test46/run.csh
#setenv PHYSLIST QGSP_BERT
=======
>>>>>>> 1.10
=======
source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QGSP_BERT_EMLSN
source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QBBC_XGGSN
source $G4INSTALL/tests/test46/run.csh
#setenv PHYSLIST QGSP_BERT
>>>>>>> 1.13
#source $G4INSTALL/tests/test46/run.csh

exit

cd $VFEM/test46
mkdir -p $REFERENCE
cd $REFERENCE

setenv PRIMARYBEAM neutron
mkdir -p $PRIMARYBEAM
cd $PRIMARYBEAM

setenv PHYSLIST QGSP_BERT_EMV
source $G4INSTALL/tests/test46/run.csh

cd $VFEM/test46
mkdir -p $REFERENCE
cd $REFERENCE

setenv PRIMARYBEAM kaon0L
mkdir -p $PRIMARYBEAM
cd $PRIMARYBEAM

setenv PRIMARYBEAM kaon0L
setenv PHYSLIST QGSP_BERT_EMV
source $G4INSTALL/tests/test46/run.csh

#

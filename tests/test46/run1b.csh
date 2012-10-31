#/bin/csh

cd $VFEM/test46
mkdir -p $REFERENCE
cd $REFERENCE

mkdir -p run1b
cd run1b

setenv PRIMARYBEAM pi-
mkdir -p $PRIMARYBEAM
cd $PRIMARYBEAM
 
setenv PHYSLIST QGSP_FTFP_BERT_EMV
#source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QGSP_BERT_EMV
#source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QBBC
source $G4INSTALL/tests/test46/run.csh

#

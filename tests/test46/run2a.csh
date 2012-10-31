#/bin/csh

cd $VFEM/test46
mkdir -p $REFERENCE
cd $REFERENCE

mkdir -p run2a
cd run2a

setenv PRIMARYBEAM proton
mkdir -p $PRIMARYBEAM
cd $PRIMARYBEAM

setenv PHYSLIST QGSP_FTFP_BERT_EMV
#source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST FTFP_BERT_EMV
source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QBBC
#source $G4INSTALL/tests/test46/run.csh

#

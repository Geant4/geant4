#/bin/csh

cd $VFEM/test46
mkdir -p $REFERENCE
cd $REFERENCE

mkdir -p run2
cd run2

setenv PRIMARYBEAM proton
mkdir -p $PRIMARYBEAM
cd $PRIMARYBEAM

setenv PHYSLIST QGSP_FTFP_BERT_EML
source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QGSP_FTFP_BERT_EML95msc93
source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QBBC
#source $G4INSTALL/tests/test46/run.csh

#

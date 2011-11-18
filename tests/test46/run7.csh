#/bin/csh

cd $VFEM/test46
mkdir -p $REFERENCE
cd $REFERENCE

setenv PRIMARYBEAM e-

mkdir -p $PRIMARYBEAM
cd $PRIMARYBEAM

setenv PHYSLIST FTFP_BERT_EMV
source $G4INSTALL/tests/test46/run_em.csh e30gev
source $G4INSTALL/tests/test46/run_em.csh e50gev

#setenv PHYSLIST QGSP_BERT_EML
#source $G4INSTALL/tests/test46/run_em.csh e50gev
#source $G4INSTALL/tests/test46/run_em.csh e30gev

#setenv PHYSLIST QGSP_BERT
#source $G4INSTALL/tests/test46/run_em.csh e50gev
#source $G4INSTALL/tests/test46/run_em.csh e30gev


#

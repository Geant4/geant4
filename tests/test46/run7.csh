#/bin/csh

cd $VFEM/test46
mkdir -p $REFERENCE

cd $REFERENCE

setenv PRIMARYBEAM e-
mkdir -p $PRIMARYBEAM
setenv PHYSLIST FTFP_BERT_EMV

cd $PRIMARYBEAM
source $G4INSTALL/tests/test46/run_em.csh e30gev
source $G4INSTALL/tests/test46/run_em.csh e50gev
source $G4INSTALL/tests/test46/run_em.csh e30gev_b
source $G4INSTALL/tests/test46/run_em.csh e50gev_b
cd ../

setenv PRIMARYBEAM gamma
mkdir -p $PRIMARYBEAM
setenv PHYSLIST FTFP_BERT_EMV

cd $PRIMARYBEAM
#source $G4INSTALL/tests/test46/run_em.csh e30gev
#source $G4INSTALL/tests/test46/run_em.csh e50gev
cd ../

#

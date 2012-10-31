#/bin/csh

cd $VFEM/test46
mkdir -p $REFERENCE

cd $REFERENCE

setenv PRIMARYBEAM pi-
mkdir -p $PRIMARYBEAM
setenv PHYSLIST QGSP_FTFP_BERT_EMV

cd $PRIMARYBEAM
source $G4INSTALL/tests/test46/run_em.csh pi-30gev
source $G4INSTALL/tests/test46/run_em.csh pi-50gev
source $G4INSTALL/tests/test46/run_em.csh pi-30gev_b
source $G4INSTALL/tests/test46/run_em.csh pi-50gev_b
cd ../

setenv PRIMARYBEAM proton
mkdir -p $PRIMARYBEAM
setenv PHYSLIST QGSP_FTFP_BERT_EMV

cd $PRIMARYBEAM
source $G4INSTALL/tests/test46/run_em.csh pi-30gev
source $G4INSTALL/tests/test46/run_em.csh pi-50gev
source $G4INSTALL/tests/test46/run_em.csh pi-30gev_b
source $G4INSTALL/tests/test46/run_em.csh pi-50gev_b
cd ../

#

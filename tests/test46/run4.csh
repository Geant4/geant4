#/bin/csh

cd $VFEM/test46
mkdir -p $REFERENCE
cd $REFERENCE

setenv PRIMARYBEAM pi-

setenv PHYSLIST QGSP_BERT
#source $G4INSTALL/tests/test46/run_em.csh e30gev100
#source $G4INSTALL/tests/test46/run_em.csh e30gev10
#source $G4INSTALL/tests/test46/run_em.csh e30gev7
#source $G4INSTALL/tests/test46/run_em.csh e30gev3
#source $G4INSTALL/tests/test46/run_em.csh e30gev1
#source $G4INSTALL/tests/test46/run_em.csh e30gev01
#source $G4INSTALL/tests/test46/run_em.csh e30gev03
#source $G4INSTALL/tests/test46/run_em.csh e30gev07
#source $G4INSTALL/tests/test46/run_em.csh e30gev001
#source $G4INSTALL/tests/test46/run_em.csh e30gev0001

setenv PHYSLIST QGSP_BERT_EMV
#source $G4INSTALL/tests/test46/run_em.csh e30gev100
#source $G4INSTALL/tests/test46/run_em.csh e30gev10
#source $G4INSTALL/tests/test46/run_em.csh e30gev7
#source $G4INSTALL/tests/test46/run_em.csh e30gev3
#source $G4INSTALL/tests/test46/run_em.csh e30gev1
#source $G4INSTALL/tests/test46/run_em.csh e30gev01
#source $G4INSTALL/tests/test46/run_em.csh e30gev03
source $G4INSTALL/tests/test46/run_em.csh e30gev07
source $G4INSTALL/tests/test46/run_em.csh e30gev001
#source $G4INSTALL/tests/test46/run_em.csh e30gev0001

#

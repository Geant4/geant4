#/bin/csh

cd $VFEM/test46
mkdir -p $REFERENCE
cd $REFERENCE

setenv PRIMARYBEAM gamma

mkdir -p $PRIMARYBEAM
cd $PRIMARYBEAM

setenv PHYSLIST QGSP_BERT
#source $G4INSTALL/tests/test46/run_em.csh g10gev100
#source $G4INSTALL/tests/test46/run_em.csh g10gev10
#source $G4INSTALL/tests/test46/run_em.csh g10gev7
#source $G4INSTALL/tests/test46/run_em.csh g10gev3
#source $G4INSTALL/tests/test46/run_em.csh g10gev1
#source $G4INSTALL/tests/test46/run_em.csh g10gev01
#source $G4INSTALL/tests/test46/run_em.csh g10gev03
#source $G4INSTALL/tests/test46/run_em.csh g10gev07
#source $G4INSTALL/tests/test46/run_em.csh g10gev001
#source $G4INSTALL/tests/test46/run_em.csh g10gev0001

setenv PHYSLIST QGSP_BERT_EML
source $G4INSTALL/tests/test46/run_em.csh g10gev100
source $G4INSTALL/tests/test46/run_em.csh g10gev10
source $G4INSTALL/tests/test46/run_em.csh g10gev7
source $G4INSTALL/tests/test46/run_em.csh g10gev3
source $G4INSTALL/tests/test46/run_em.csh g10gev1
source $G4INSTALL/tests/test46/run_em.csh g10gev01
source $G4INSTALL/tests/test46/run_em.csh g10gev03
source $G4INSTALL/tests/test46/run_em.csh g10gev07
source $G4INSTALL/tests/test46/run_em.csh g10gev001
source $G4INSTALL/tests/test46/run_em.csh g10gev0001

#

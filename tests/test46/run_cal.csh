#/bin/csh

cd $VFEM/test46
mkdir -p $REFERENCE

cd $REFERENCE

setenv PRIMARYBEAM pi-
mkdir -p $PRIMARYBEAM
cd $PRIMARYBEAM

setenv PHYSLIST QGSP_FTFP_BERT_EML
source $G4INSTALL/tests/test46/run_em.csh e50gev_cal_b

setenv PHYSLIST QGSP_FTFP_BERT_EML95msc93
source $G4INSTALL/tests/test46/run_em.csh e50gev_cal_b

setenv PHYSLIST FTFP_BERT_EMV
source $G4INSTALL/tests/test46/run_em.csh e50gev_cal_b


#

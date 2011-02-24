#/bin/csh

cd $VFEM/test46
mkdir -p $REFERENCE
cd $REFERENCE

setenv PRIMARYBEAM proton
mkdir -p $PRIMARYBEAM
cd $PRIMARYBEAM

setenv PHYSLIST QBBC
#source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QBBCG
#source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QBBCF
#source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QGSP_BERT_EML
<<<<<<< run2.csh
<<<<<<< run2.csh
source $G4INSTALL/tests/test46/run.csh
setenv PHYSLIST QGSP_BERT_EMLSN
source $G4INSTALL/tests/test46/run.csh
=======
#source $G4INSTALL/tests/test46/run.csh
=======
source $G4INSTALL/tests/test46/run.csh
>>>>>>> 1.13
setenv PHYSLIST QGSP_BERT_EMV
<<<<<<< run2.csh
#source $G4INSTALL/tests/test46/run.csh
>>>>>>> 1.10
=======
source $G4INSTALL/tests/test46/run.csh
>>>>>>> 1.13
setenv PHYSLIST FTFP_BERT
#source $G4INSTALL/tests/test46/run.cs
setenv PHYSLIST QGSP_BERT
#source $G4INSTALL/tests/test46/run.csh

#

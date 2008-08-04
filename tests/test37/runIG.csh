#/bin/csh

mkdir -p $REFERENCE
cd $REFERENCE
rm *

set    work = "$G4MY/test37"
set    dir  = "$G4INSTALL/tests/test37/"

setenv PHYSLIST    standardIG
set    phys = "optG"
source ${dir}run_single.csh ${phys} ${work} ${dir}

ln -s $VFEM/test37/geant4-09-01-ref-03/*S.* ./
ln -s $VFEM/test37/geant4-09-01-ref-03/*0.* ./
ln -s $VFEM/test37/geant4-09-01-ref-03/*1.* ./
#ln -s $VFEM/test37/geant4-09-01-ref-03/*3.* ./
ln -s ../20080804/*3.* ./

$G4MY/reader_test37 Al     $1
$G4MY/reader_test37 Mo     $1
$G4MY/reader_test37 Ta     $1
$G4MY/reader_test37 TaAl   $1
$G4MY/reader_test37 AlAuAl $1

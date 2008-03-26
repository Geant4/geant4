#/bin/csh

mkdir -p $REFERENCE
cd $REFERENCE

set    work = "$G4MY/test37"
set    dir  = "$G4INSTALL/tests/test37/"

setenv PHYSLIST    emstandard
set    phys = "opt0"
source ${dir}run_single.csh ${phys} ${work} ${dir}

setenv PHYSLIST    standard
set    phys = "opt3"
source ${dir}run_single.csh ${phys} ${work} ${dir}

setenv PHYSLIST    emstandard_opt1
set    phys = "opt1"
source ${dir}run_single.csh ${phys} ${work} ${dir}

setenv PHYSLIST    standardIG
set    phys = "optG"
source ${dir}run_single.csh ${phys} ${work} ${dir}

cp $VFEM/test37/geant4-V09-01-ref-00/*S.log ./

$G4MY/reader_test37 Al     $1
$G4MY/reader_test37 Mo     $1
$G4MY/reader_test37 Ta     $1
$G4MY/reader_test37 TaAl   $1
$G4MY/reader_test37 AlAuAl $1

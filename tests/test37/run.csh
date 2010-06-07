#/bin/csh

cd $VFEM/test37
mkdir -p $REFERENCE
cd $REFERENCE
rm *

set    work = "$G4MY/test37"
set    dir  = "$G4INSTALL/tests/test37/"

setenv PHYSLIST    emstandard
set    phys = "opt0"
source ${dir}run_single.csh ${phys} ${work} ${dir}

setenv PHYSLIST    emstandard_opt3
set    phys = "opt3"
source ${dir}run_single.csh ${phys} ${work} ${dir}

setenv PHYSLIST    standardGS
set    phys = "optK"
source ${dir}run_single.csh ${phys} ${work} ${dir}

setenv PHYSLIST    standard_opt2
set    phys = "optG"
source ${dir}run_single.csh ${phys} ${work} ${dir}

setenv PHYSLIST    standardSS
set    phys = "optS"
source ${dir}run_single.csh ${phys} ${work} ${dir}

$G4MY/reader_test37 Al     $1
$G4MY/reader_test37 Mo     $1
$G4MY/reader_test37 Ta     $1
$G4MY/reader_test37 TaAl   $1
$G4MY/reader_test37 AlAuAl $1

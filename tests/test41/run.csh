#/bin/csh

mkdir -p $REFERENCE
cd $REFERENCE

set    work = "$G4MY/test41"
set    dir  = "$G4INSTALL/tests/test41/"

setenv PHYSLIST    emstandard
set    phys = "opt0"
source ${dir}run_single.csh ${phys} ${work} ${dir}

setenv PHYSLIST    standard
set    phys = "opt3"
source ${dir}run_single.csh ${phys} ${work} ${dir}

setenv PHYSLIST    emstandard_opt1
set    phys = "opt1"
source ${dir}run_single.csh ${phys} ${work} ${dir}

setenv PHYSLIST    standardSS
set    phys = "optS"
source ${dir}run_single.csh ${phys} ${work} ${dir}

source plot.csh $1

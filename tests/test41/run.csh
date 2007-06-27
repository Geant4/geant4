#/bin/csh

mkdir -p $REFERENCE
cd $REFERENCE

set    work = "$G4MY/test41 $G4INSTALL/tests/test41/"

setenv PHYSLIST    emstandard
set    phys = "opt0"
source run_single.csh ${phys} ${work}

setenv PHYSLIST    standard
set    phys = "opt3"
source run_single.csh ${phys} ${work}

setenv PHYSLIST    standard_opt1
set    phys = "opt1"
source run_single.csh ${phys} ${work}

setenv PHYSLIST    standardSS
set    phys = "optS"
source run_single.csh ${phys} ${work}

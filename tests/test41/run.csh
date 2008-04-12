#/bin/csh

mkdir -p $REFERENCE
cd $REFERENCE
rm -f r.out p.out

set    work = "$G4MY/test41"
set    dir  = "$G4INSTALL/tests/test41/"

setenv PHYSLIST    emstandard
set    phys = "opt0"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& res0.out

setenv PHYSLIST    standard_local
set    phys = "opt3"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& res3.out

setenv PHYSLIST    emstandard_opt1
set    phys = "opt1"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& res1.out

setenv PHYSLIST    standardIG
set    phys = "optG"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& resG.out

setenv PHYSLIST    standardSS
set    phys = "optS"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& resS.out

source ${dir}plot.csh $1 >& p.out

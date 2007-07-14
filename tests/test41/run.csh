#/bin/csh

mkdir -p $REFERENCE
cd $REFERENCE
rm -f r.out p.out

set    work = "$G4MY/test41"
set    dir  = "$G4INSTALL/tests/test41/"

setenv PHYSLIST    emstandard
set    phys = "opt0"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& r.out

setenv PHYSLIST    standard
set    phys = "opt3"
source ${dir}run_single.csh ${phys} ${work} ${dir} >>& r.out

setenv PHYSLIST    emstandard_opt1
set    phys = "opt1"
source ${dir}run_single.csh ${phys} ${work} ${dir} >>& r.out

setenv PHYSLIST    standardSS
set    phys = "optS"
source ${dir}run_single.csh ${phys} ${work} ${dir} >>& r.out

source ${dir}plot.csh $1 >& p.out

#/bin/csh

mkdir -p $REFERENCE
cd $REFERENCE
rm -f r.out p.out

set    work = "$G4MY/test41"
set    dir  = "$G4INSTALL/tests/test41/"

setenv PHYSLIST    emstandard
set    phys = "opt0"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& r.out

setenv PHYSLIST    standard_local
set    phys = "opt3"
source ${dir}run_single.csh ${phys} ${work} ${dir} >>& r.out

setenv PHYSLIST    emstandard_opt1
set    phys = "opt1"
source ${dir}run_single.csh ${phys} ${work} ${dir} >>& r.out

setenv PHYSLIST    standardIG
set    phys = "optG"
source ${dir}run_single.csh ${phys} ${work} ${dir} >>& r.out

cp $VFEM/test41/geant4-09-01-ref-03/*S.log ./
source ${dir}plot.csh $1 >& p.out

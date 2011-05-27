#/bin/csh

cd $VFEM/test41
mkdir -p $REFERENCE
cd $REFERENCE
rm -f *.*

set    work = "$G4MY/test41"
set    dir  = "$G4INSTALL/tests/test41/"

setenv PHYSLIST    emstandard
set    phys = "opt0"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& res0.out

setenv PHYSLIST    emstandard_msc95
set    phys = "msc95"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& res2.out

setenv PHYSLIST    emstandard_msc93
set    phys = "msc93"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& res93.out

setenv PHYSLIST    emstandard_msc90
set    phys = "msc90"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& res3.out

setenv PHYSLIST    standardSS
set    phys = "optSS"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& resS.out

source ${dir}plot.csh $1 >& p.out

gzip *.out

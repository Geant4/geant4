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
source ${dir}run_single.csh ${phys} ${work} ${dir} >& res95.out

setenv PHYSLIST    emstandard_msc93
set    phys = "msc93"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& res93.out

setenv PHYSLIST    emstandard_msc96
set    phys = "msc96"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& res96.out

setenv PHYSLIST    standardSS
set    phys = "optSS"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& resSS.out

source ${dir}plot.csh $1 >& p.out
$G4MY/reader_test41 $1 >& p1.out

gzip *.out

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

setenv PHYSLIST    emstandard_opt3
set    phys = "opt3"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& res3.out

setenv PHYSLIST    emstandardIG
set    phys = "mscWVI"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& resWVI.out

setenv PHYSLIST    emstandardGS
set    phys = "mscGS"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& resGS.out

setenv PHYSLIST    standardSS
set    phys = "optSS"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& resSS.out

source ${dir}plot.csh $1 >& p.out
$G4MY/reader_test41 $1 >& p1.out

gzip *.out

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

setenv PHYSLIST    emstandardUB
set    phys = "optUB"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& resUB.out

setenv PHYSLIST    emstandardIG
set    phys = "mscWVI"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& resWVI.out

setenv PHYSLIST    emstandardSS
set    phys = "optSS"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& resSS.out

setenv PHYSLIST    emstandardSSM
set    phys = "optSSM"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& resSSM.out

source ${dir}plot.csh $1 >& p.out
$G4MY/reader_test41 $1 >& p1.out

gzip *.out

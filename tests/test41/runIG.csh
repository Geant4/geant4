#/bin/csh

cd $VFEM/test41
mkdir -p $REFERENCE
cd $REFERENCE
rm *

set    work = "$G4MY/test41"
set    dir  = "$G4INSTALL/tests/test41/"

setenv PHYSLIST    emstandard
set    phys = "opt0"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& res0.out

ln -s $VFEM/test41/$2/*S.* ./
ln -s $VFEM/test41/$2/*96.* ./
ln -s $VFEM/test41/$2/*93.* ./
ln -s $VFEM/test41/$2/*95.* ./

source ${dir}plot.csh $1 >& p.out
$G4MY/reader_test41 $1 >& p1.out
#

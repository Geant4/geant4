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

ln -s $VFEM/test41/$2/*SS.* ./
ln -s $VFEM/test41/$2/*SSM.* ./
ln -s $VFEM/test41/$2/*UB.* ./
ln -s $VFEM/test41/$2/*WVI.* ./

source ${dir}plot.csh $1 >& p.out
$G4MY/reader_test41 $1 >& p1.out
#

#/bin/csh

cd $VFEM/test41
mkdir -p $REFERENCE
cd $REFERENCE
rm *

set    work = "$G4MY/test41"
set    dir  = "$G4INSTALL/tests/test41/"

setenv PHYSLIST    standardIG
set    phys = "opt2"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& resG.out

ln -s $VFEM/test41/$2/*S.* ./
ln -s $VFEM/test41/$2/*0.* ./
ln -s $VFEM/test41/$2/*3.* ./

source ${dir}plot.csh $1 >& p.out

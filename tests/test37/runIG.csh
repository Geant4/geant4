#/bin/csh

cd $VFEM/test37
mkdir -p $REFERENCE
cd $REFERENCE
rm *

set    work = "$G4MY/test37"
set    dir  = "$G4INSTALL/tests/test37/"

setenv PHYSLIST    emstandardWVI
set    phys = "optG"
source ${dir}run_single.csh ${phys} ${work} ${dir}

ln -s $VFEM/test37/$2/*S.* ./
ln -s $VFEM/test37/$2/*0.* ./
ln -s $VFEM/test37/$2/*3.* ./
ln -s $VFEM/test37/$2/*K.* ./

$G4MY/reader_test37 Al     $1
$G4MY/reader_test37 Be     $1
$G4MY/reader_test37 Mo     $1
$G4MY/reader_test37 Ta     $1
$G4MY/reader_test37 TaAl   $1
$G4MY/reader_test37 AlAuAl $1

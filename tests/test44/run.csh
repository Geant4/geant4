#/bin/csh

cd $VFEM/test44

mkdir -p $REFERENCE
cd $REFERENCE
rm -f *.txt

set    work = "$G4BIN/$G4SYSTEM/test44"
set    dir  = "$G4INSTALL/tests/test44/"

ln -s ${dir}/Exp_Data/*.txt ./

setenv PHYSLIST  QBBC
set    phys = "opt0"
source ${dir}run_single.csh ${phys} ${work} ${dir}

set    phys = "opt2"
source ${dir}run_single.csh ${phys} ${work} ${dir}

set    phys = "opt3"
source ${dir}run_single.csh ${phys} ${work} ${dir}

${G4BIN}/${G4SYSTEM}/reader_test44 p $1
${G4BIN}/${G4SYSTEM}/reader_test44 he4 $1
${G4BIN}/${G4SYSTEM}/reader_test44 c12 $1

source ${dir}plot.csh $1 >& root.log


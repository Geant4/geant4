#/bin/csh

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

source plot.csh $1 >& root.log

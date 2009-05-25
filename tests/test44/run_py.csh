#/bin/csh

if ( $?LCG_SYS == 0 ) then
setenv LCG_SYS    slc4_ia32_gcc34
endif

if ($?LD_LIBRARY_PATH == 0) then
setenv LD_LIBRARY_PATH /usr/local/lib
endif

if ( $?PYTHONPATH == 0 ) then
setenv PYTHONPATH /afs/cern.ch/sw/lcg/external/Python/2.5.4/${LCG_SYS}
endif

if ( $?ROOTSYS == 0 ) then
setenv ROOTSYS /afs/cern.ch/sw/lcg/app/releases/ROOT/5.21.06/$LCG_SYS/root
endif

setenv PYTHONDIR /afs/cern.ch/sw/lcg/external/Python/2.5.4/${LCG_SYS}
setenv PATH $PYTHONDIR/bin:$PATH
setenv LD_LIBRARY_PATH $PYTHONDIR/lib:$LD_LIBRARY_PATH
setenv PYTHONPATH $ROOTSYS/lib:$PYTHONPATH


if ( $?REFERENCE == 0 ) then
setenv REFERENCE `date '+%m_%d_%Y-%H:%M:%S'`
endif

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

chmod +x ../utils/reader_test44.py
../utils/reader_test44.py p $1
../utils/reader_test44.py he4 $1
../utils/reader_test44.py c12 $1

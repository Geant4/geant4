#/bin/csh
#trace on
#---------------------------------------------------------------
#
#  Author V.Ivanchenko 27 May 2008
#  test45
#
#----------------------------------------------------------------


mkdir -p $REFERENCE
cd $REFERENCE
rm -f *.TXT

set    work = "$G4BIN/$G4SYSTEM/test44"
set    dir  = "$G4INSTALL/tests/test44/"

ln -f ${dir}/Exp_Data/*.txt

setenv PHYSLIST    QBBC
set    phys = "var0"
source ${dir}run_single.csh ${phys} ${work} ${dir}

#source ${dir}plot.csh $1
#source ${dir}plot.csh 
#${dir}/utils/reader_test44  p_water_${phys}.log
${G4BIN}/${G4SYSTEM}/reader_test44 ${phys} $1

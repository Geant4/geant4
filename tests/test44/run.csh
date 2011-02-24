#/bin/csh

if ( $?REFERENCE == 0 ) then
setenv REFERENCE `date '+%m_%d_%Y-%H:%M:%S'`
endif

cd $VFEM/test44

mkdir -p $REFERENCE
cd $REFERENCE
rm -f *.txt *.gz

set    work = "$G4BIN/$G4SYSTEM/test44"
set    dir  = "$G4INSTALL/tests/test44/"

ln -s ${dir}/Exp_Data/*.txt ./

setenv PHYSLIST  QBBC

set tPart = (p he4 c12)
foreach phys (opt0 opt2 opt3) 
    foreach part ($tPart)
	set file = "${tPart}_water_${phys}.log"
	if ( -e "$file" )  then
	    rm -f $file
	endif
	${work} ${dir}${part}_water_${phys}.in >& ${part}_water_${phys}.log
	mv Bragg.out ${part}_${phys}.out
    end
end

date > root.log
foreach tPart ($tPart)
    ${G4BIN}/${G4SYSTEM}/reader_test44 $tPart $1  >>& root.log
end

gzip *.log *.out

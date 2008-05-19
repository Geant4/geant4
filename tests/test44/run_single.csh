#/bin/csh

rm -f p_water_${1}.log
${2} ${3}p_water.in >& p_water_${1}.log
mv Bragg.out ${1}.out

#




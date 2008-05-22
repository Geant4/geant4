#/bin/csh

rm -f p_water_${1}.log he4_water_${1}.log c12_water_${1}.log

${2} ${3}p_water_${1}.in >& p_water_${1}.log
mv Bragg.out p_${1}.out

${2} ${3}he4_water_${1}.in >& he4_water_${1}.log
mv Bragg.out he4_${1}.out

${2} ${3}c12_water_${1}.in >& c12_water_${1}.log
mv Bragg.out c12_${1}.out

#




#/bin/csh

rm -f ${1}.out

echo 'Run test37 at ' $HOST >& $1.out

${2} ${3}Aluminum.in >>& ${1}.out
mv Sandia.out Al_${1}.log

${2} ${3}Molybdenum.in >>& ${1}.out
mv Sandia.out Mo_${1}.log

${2} ${3}Tantalum.in >>& ${1}.out
mv Sandia.out Ta_${1}.log

${2} ${3}TaAl.in >>& ${1}.out
mv Sandia.out TaAl_${1}.log

${2} ${3}AlAuAl.in >>& ${1}.out
mv Sandia.out AlAuAl_${1}.log

${2} ${3}Beryllium.in >>& ${1}.out
mv Sandia.out Be_${1}.log

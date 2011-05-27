#!/bin/csh -f
# Main file for production 

set i = 1
set F = ( QGSP_BERT FTFP_BERT LHEP CHIPS QGSP_BIC QGSP_INCL_ABLA Shielding )

mkdir -p $RESULTS/HADR02/$REFERENCE/ 


while ($i < 8)
setenv PHYSLIST $F[$i]
mkdir -p $RESULTS/HADR02/$REFERENCE/$F[$i]/
cd $RESULTS/HADR02/$REFERENCE/$F[$i]/
echo "PRODUCTION in " $RESULTS/HADR02/$REFERENCE/$F[$i] 
echo $F[$i]

$G4WORKDIR/bin/Linux-g++/hadr02 $HADR02/dpmjet1.in >&! r.out 

@ i++
echo " "
end 
exit






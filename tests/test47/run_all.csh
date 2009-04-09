#!/bin/tcsh -f

# this is CMSLPC setup
#
###cd /uscms_data/d2/yarba_j/geant4/tests
#
# this is LXPLUS
#
cd $HOME/scratch0/geant4/tests/test47

./run_itep.csh 12 0

set i = 0
set njobs = 3
while ( $i < $njobs ) 
/usr/bin/bsub -q 2nd run_bnl.csh $i
@ i = $i + 1
end


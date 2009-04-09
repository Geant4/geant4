#!/bin/tcsh -f

# this is CMSLPC setup
#
###cd /uscms_data/d2/yarba_j/geant4/tests
#
# this is LXPLUS
#
cd $HOME/scratch0/geant4/tests/test47


if ( $#argv == 0 ) then
set i = 0
set njobs = 12
else
   if ( $#argv == 1 ) then
      set i = 0
      set njobs = $1
   else
      set i = $2
      set njobs = $1
   endif
endif


while ( $i < $njobs )

if ( $i < 4 ) then
/usr/bin/bsub -q 2nd run_pip_itep.csh $i
else
  if ( $i < 8 ) then
    /usr/bin/bsub -q 2nd run_pim_itep.csh $i
  else
    /usr/bin/bsub -q 2nd run_proton_itep.csh $i
  endif
endif 

@ i = $i + 1

end

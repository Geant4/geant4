#!/bin/csh -f
#trace on
#---------------------------------------------------------------
#
#  Author V.Ivanchenko 27 April 2004
#  Standard EM test for pcgeant4
#
#----------------------------------------------------------------
  
  cd ${G4INSTALL}/source
  gmake
  cd ../examples/extended/electromagnetic

  set ncurr  = 0
    while ($ncurr <= 14)
      cd TestEm$ncurr
      gmake
      $G4MY/TestEm$ncurr TestEm$ncurr.in >& $1
      cd ../
      @ ncurr++
    end

  set ncurr  = 0
    while ($ncurr <= 14)
      tkdiff TestEm$ncurr/$2  TestEm$ncurr/$1
      @ ncurr++
    end

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

  set ncurr  = 1
    while ($ncurr <= 10)
      cd TestEm$ncurr
      gmake
      $G4MY/TestEm$ncurr TestEm$ncurr.in >& $1
      cd ../
      @ ncurr++
    end
    cd PhotonProcesses
    gmake 
    $G4MY/PhotonProcesses PhotonProcesses.in >& $1
    cd ../MuonProcesses
    gmake 
    $G4MY/MuonProcesses MuonProcesses.in >& $1
    cd ../

  set ncurr  = 1
    while ($ncurr <= 10)
      tkdiff TestEm$ncurr/$2  TestEm$ncurr/$1
      @ ncurr++
    end
    tkdiff PhotonProcesses/$2  PhotonProcesses/$1
    tkdiff MuonProcesses/$2  MuonProcesses/$1

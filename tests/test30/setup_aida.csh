ifndef G4ANALYSIS_BUILD
  setenv G4ANALYSIS_BUILD  1
  source /afs/cern.ch/sw/lhcxx/share/LHCXX/4.0.4/install/sharedstart.csh

  setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${OGLHOME}/lib 
endif

echo $LD_LIBRARY_PATH

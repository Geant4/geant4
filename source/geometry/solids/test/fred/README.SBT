$Id: README.SBT,v 1.1 2008-02-22 11:20:39 ivana Exp $

The SBT log outputs can be debugged with Fred in this
way:

1. Copy a selected SBT log output file in Fred/test3.log

   $>  cd $G4INSTALL/source/geometry/solids/test/fred
   $>  cp ../SBT/log/cons.e2.log test3.log

2. Update the solid definition in src/FredDetectorConstruction.cc
   and recompile Fred:
   
   $> make

3, Uncomment the solid type in test1.mac   
   /fred/volume CONE

4. Run fred and process all errors in the log file:

   $> fred
   /control/execute test1.mac
   /fred/test3/runDebug

   Or debug a selected error using the commands: 

   /fred/test3/debug 3
   /run/beamOn 1
  
        

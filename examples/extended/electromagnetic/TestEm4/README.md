\page ExampleTestEm4 Example TestEm4

 Plot energy deposited by 9 MeV photon beam in an homogeneous medium.
	
## GEOMETRY DEFINITION
 
 It is a cylinder of 5 cm radius filled with C6F6.
 	
## PHYSICS LIST
 
 The particle list contains only gamma, electron,positron.
 The physics list contains the 'standard' electromagnetic processes.
 	 
## AN EVENT : THE PRIMARY GENERATOR
 
 The primary kinematic is a single 9 MeV gamma randomly shooted at the
 middle of the cylinder. 
 	 				
## VISUALIZATION
 
 The Visualization Manager is set in the main().
 The initialisation of the drawing is done via the commands
 `/vis/..` in the macro vis.mac. This macro is
 automatically read from the main in case of interactive running mode.
 
 The detector has a default view which is a transversal view of the 
 cylinder.
 
 The tracks are drawn at the end of event, and erased at the end of run.
 Optionaly one can choose to draw all particles, only the charged one,
 or none. This command is defined in EventActionMessenger class.
 	
## PHYSICS SURVEY
 
 The energy deposited in C6F6 is histogramed.
 	
## HOW TO START ?
 
 - Execute TestEm4 in 'batch' mode from macro files
```
% ./TestEm4   TestEm4.in
```
  		
 - Execute TestEm4 in 'interactive mode' with visualization
```
% ./TestEm4
....
Idle> type your commands
....
Idle> exit
```

macro verbose.mac illustrate capability of verbosity
  		
## USING HISTOGRAMS

 The format of the histogram file can be : root (default),
 xml, csv, by changing the default file type in RunAction.cc

## RANDOM NUMBERS HANDLING
 
   CLHEP provides several random number engines. In this example the Ranecu
   engine is choosen at beginning of the main (TestEm4.cc).
	
   By default, G4RunManager does not save the rndm seed.
   To do so the user must set in BeginOfRunAction:
   G4RunManager::GetRunManager()->SetRandomNumberStore(true);
	
   Then the rndm seed is systematically saved at beginning of run
   (currentRun.rndm) and beginning of event (currentEvent.rndm)
   Therefore, in case of abnormal end, the seed of the last event processed
   is available in currentEvent.rndm
	
   Even in case of normal run processing, the user may wish to preserve the
   rndm seed of selected events. At any time in the event, put the
   following statement:
```cpp
if (condition) G4RunManager::GetRunManager()->rndmSaveThisEvent();
```
   currentEvent.rndm will be copied to runXXevntYY.rndm
   (see SteppingAction::UserSteppingAction() )
	
   To restart a run from a given rndm seed, use the UI command :
```
/random/resetEngineFrom  fileName 
```
	
   The macro rndmSeed.mac shows how to save and reset the random number
   seed between runs, from UI commands.	

  
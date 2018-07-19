//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4RunMessenger.cc 95634 2016-02-17 08:05:21Z gcosmo $
//

#include "G4RunMessenger.hh"
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UImanager.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ios.hh"
#include "G4MaterialScanner.hh"
#include "G4Tokenizer.hh"
#include "Randomize.hh"
#include <sstream>

G4RunMessenger::G4RunMessenger(G4RunManager * runMgr)
:runManager(runMgr)
{
  runDirectory = new G4UIdirectory("/run/");
  runDirectory->SetGuidance("Run control commands.");

  initCmd = new G4UIcmdWithoutParameter("/run/initialize",this);
  initCmd->SetGuidance("Initialize G4 kernel.");
  initCmd->AvailableForStates(G4State_PreInit,G4State_Idle),

  beamOnCmd = new G4UIcommand("/run/beamOn",this);
  beamOnCmd->SetGuidance("Start a Run.");
  beamOnCmd->SetGuidance("If G4 kernel is not initialized, it will be initialized.");
  beamOnCmd->SetGuidance("Default number of events to be processed is 1.");
  beamOnCmd->SetGuidance("The second and third arguments can be used for");
  beamOnCmd->SetGuidance("executing a macro file at the end of each event.");
  beamOnCmd->SetGuidance("If the second argument, i.e. name of the macro");
  beamOnCmd->SetGuidance("file, is given but the third argument is not,");
  beamOnCmd->SetGuidance("the macro file will be executed for all of the");
  beamOnCmd->SetGuidance("event.");
  beamOnCmd->SetGuidance("If the third argument (nSelect) is given, the");
  beamOnCmd->SetGuidance("macro file will be executed only for the first");
  beamOnCmd->SetGuidance("nSelect events.");
  beamOnCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  beamOnCmd->SetToBeBroadcasted(false);
  G4UIparameter* p1 = new G4UIparameter("numberOfEvent",'i',true);
  p1->SetDefaultValue(1);
  p1->SetParameterRange("numberOfEvent >= 0");
  beamOnCmd->SetParameter(p1);
  G4UIparameter* p2 = new G4UIparameter("macroFile",'s',true);
  p2->SetDefaultValue("***NULL***");
  beamOnCmd->SetParameter(p2);
  G4UIparameter* p3 = new G4UIparameter("nSelect",'i',true);
  p3->SetDefaultValue(-1);
  p3->SetParameterRange("nSelect>=-1");
  beamOnCmd->SetParameter(p3);
  //beamOnCmd->SetToBeBroadcasted(false);

  verboseCmd = new G4UIcmdWithAnInteger("/run/verbose",this);
  verboseCmd->SetGuidance("Set the Verbose level of G4RunManager.");
  verboseCmd->SetGuidance(" 0 : Silent (default)");
  verboseCmd->SetGuidance(" 1 : Display main topics");
  verboseCmd->SetGuidance(" 2 : Display main topics and run summary");
  verboseCmd->SetParameterName("level",true);
  verboseCmd->SetDefaultValue(0);
  verboseCmd->SetRange("level >=0 && level <=2");

  printProgCmd = new G4UIcmdWithAnInteger("/run/printProgress",this);
  printProgCmd->SetGuidance("Display begin_of_event information at given frequency.");
  printProgCmd->SetGuidance("If it is set to zero, only the begin_of_run is shown.");
  printProgCmd->SetGuidance("If it is set to -1, no print-out is shown.");
  printProgCmd->SetParameterName("mod",true);
  printProgCmd->SetDefaultValue(-1);
  printProgCmd->SetRange("mod>=-1");

  nThreadsCmd = new G4UIcmdWithAnInteger("/run/numberOfThreads",this);
  nThreadsCmd->SetGuidance("Set the number of threads to be used.");
  nThreadsCmd->SetGuidance("This command works only in PreInit state.");
  nThreadsCmd->SetGuidance("This command is valid only for multi-threaded mode.");
  nThreadsCmd->SetGuidance("The command is ignored if it is issued in sequential mode.");
  nThreadsCmd->SetParameterName("nThreads",true);
  nThreadsCmd->SetDefaultValue(2);
  nThreadsCmd->SetRange("nThreads >0");
  nThreadsCmd->SetToBeBroadcasted(false);
  nThreadsCmd->AvailableForStates(G4State_PreInit);

  maxThreadsCmd = new G4UIcmdWithoutParameter("/run/useMaximumLogicalCores",this);
  maxThreadsCmd->SetGuidance("Set the number of threads to be the number of available logical cores.");
  maxThreadsCmd->SetGuidance("This command works only in PreInit state.");
  maxThreadsCmd->SetGuidance("This command is valid only for multi-threaded mode.");
  maxThreadsCmd->SetGuidance("The command is ignored if it is issued in sequential mode.");
  maxThreadsCmd->SetToBeBroadcasted(false);
  maxThreadsCmd->AvailableForStates(G4State_PreInit);

  pinAffinityCmd = new G4UIcmdWithAnInteger("/run/pinAffinity",this);
  pinAffinityCmd->SetGuidance("Locks each thread to a specific logical core. Workers are locked in round robin to logical cores.");
  pinAffinityCmd->SetGuidance("This command is valid only for multi-threaded mode.");
  pinAffinityCmd->SetGuidance("This command works only in PreInit state.");
  pinAffinityCmd->SetGuidance("This command is ignored if it is issued in sequential mode.");
  pinAffinityCmd->SetGuidance("If a value n>0 is provided it starts setting affinity from the n-th CPU (note: counting from 1).");
  pinAffinityCmd->SetGuidance("E.g. /run/pinAffinity 3 locks first thread on third logical CPU (number 2).");
  pinAffinityCmd->SetGuidance("If a value n<0 is provided never locks on n-th CPU.");
  pinAffinityCmd->SetParameterName("pinAffinity",true);
  pinAffinityCmd->SetDefaultValue(1);
  pinAffinityCmd->SetToBeBroadcasted(false);
  pinAffinityCmd->SetRange("pinAffinity > 0 || pinAffinity < 0");
  pinAffinityCmd->AvailableForStates(G4State_PreInit);

  evModCmd = new G4UIcommand("/run/eventModulo",this);
  evModCmd->SetGuidance("Set the event modulo for dispatching events to worker threads"); 
  evModCmd->SetGuidance("i.e. each worker thread is ordered to simulate N events and then");
  evModCmd->SetGuidance("comes back to G4MTRunManager for next set.");
  evModCmd->SetGuidance("If it is set to zero (default value), N is roughly given by this.");
  evModCmd->SetGuidance("   N = int( sqrt( number_of_events / number_of_threads ) )");
  evModCmd->SetGuidance("The value N may affect on the computing performance in particular");
  evModCmd->SetGuidance("if N is too small compared to the total number of events.");
  evModCmd->SetGuidance("The second parameter seedOnce specifies how frequently each worker");
  evModCmd->SetGuidance("thread is seeded by the random number sequence contrally managed");
  evModCmd->SetGuidance("by the master G4MTRunManager.");
  evModCmd->SetGuidance(" - If seedOnce is set to 0 (default), seeds that are centrally managed");
  evModCmd->SetGuidance("   by G4MTRunManager are set for every event of every worker thread.");
  evModCmd->SetGuidance("   This option guarantees event reproducability regardless of number");
  evModCmd->SetGuidance("   of threads.");
  evModCmd->SetGuidance(" - If seedOnce is set to 1, seeds are set only once for the first");
  evModCmd->SetGuidance("   event of each run of each worker thread. Event reproducability is");
  evModCmd->SetGuidance("   guaranteed only if the same number of worker threads are used.");
  evModCmd->SetGuidance("   On the other hand, this option offers better computing performance");
  evModCmd->SetGuidance("   in particular for applications with relatively small primary");
  evModCmd->SetGuidance("   particle energy and large number of events."); 
  evModCmd->SetGuidance(" - If seedOnce is set to 2, seeds are set only for the first event of");
  evModCmd->SetGuidance("   group of N events. This option is reserved for the future use when");
  evModCmd->SetGuidance("   Geant4 allows number of threads to be dynatically changed during an");
  evModCmd->SetGuidance("   event loop.");
  evModCmd->SetGuidance("This command is valid only for multi-threaded mode.");
  evModCmd->SetGuidance("This command is ignored if it is issued in sequential mode.");
  G4UIparameter* emp1 = new G4UIparameter("N",'i',true);
  emp1->SetDefaultValue(0);
  emp1->SetParameterRange("N >= 0");
  evModCmd->SetParameter(emp1);
  G4UIparameter* emp2 = new G4UIparameter("seedOnce",'i',true);
  emp2->SetDefaultValue(0);
  emp2->SetParameterRange("seedOnce >= 0 && seedOnce <=2");
  evModCmd->SetParameter(emp2);
  evModCmd->SetToBeBroadcasted(false);
  evModCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  dumpRegCmd = new G4UIcmdWithAString("/run/dumpRegion",this);
  dumpRegCmd->SetGuidance("Dump region information.");
  dumpRegCmd->SetGuidance("In case name of a region is not given, all regions will be displayed.");
  dumpRegCmd->SetParameterName("regionName", true);
  dumpRegCmd->SetDefaultValue("**ALL**");
  dumpRegCmd->AvailableForStates(G4State_Idle);

  dumpCoupleCmd = new G4UIcmdWithoutParameter("/run/dumpCouples",this);
  dumpCoupleCmd->SetGuidance("Dump material-cuts-couple information.");
  dumpCoupleCmd->SetGuidance("Note that material-cuts-couple information is updated");
  dumpCoupleCmd->SetGuidance("after BeamOn has started.");
  dumpCoupleCmd->AvailableForStates(G4State_Idle);

  optCmd = new G4UIcmdWithABool("/run/optimizeGeometry",this);
  optCmd->SetGuidance("Set the optimization flag for geometry.");
  optCmd->SetGuidance("If it is set to TRUE, G4GeometryManager will optimize");
  optCmd->SetGuidance("the geometry definitions.");
  optCmd->SetGuidance("GEANT4 is initialized with this flag as TRUE.");
  optCmd->SetParameterName("optimizeFlag",true);
  optCmd->SetDefaultValue(true);
  optCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  brkBoECmd = new G4UIcmdWithABool("/run/breakAtBeginOfEvent",this);
  brkBoECmd->SetGuidance("Set a break point at the begining of every event.");
  brkBoECmd->SetParameterName("flag",true);
  brkBoECmd->SetDefaultValue(true);
  
  brkEoECmd = new G4UIcmdWithABool("/run/breakAtEndOfEvent",this);
  brkEoECmd->SetGuidance("Set a break point at the end of every event.");
  brkEoECmd->SetParameterName("flag",true);
  brkEoECmd->SetDefaultValue(true);
  
  abortCmd = new G4UIcmdWithABool("/run/abort",this);
  abortCmd->SetGuidance("Abort current run processing.");
  abortCmd->SetGuidance("If softAbort is false (default), currently processing event will be immediately aborted,");
  abortCmd->SetGuidance("while softAbort is true, abortion occurs after processing the current event.");
  abortCmd->AvailableForStates(G4State_GeomClosed,G4State_EventProc);
  abortCmd->SetParameterName("softAbort",true);
  abortCmd->SetDefaultValue(false);

  abortEventCmd = new G4UIcmdWithoutParameter("/run/abortCurrentEvent",this);
  abortEventCmd->SetGuidance("Abort currently processing event.");
  abortEventCmd->AvailableForStates(G4State_EventProc);

  geomCmd = new G4UIcmdWithoutParameter("/run/geometryModified",this);
  geomCmd->SetGuidance("Force geometry to be closed (re-voxellized) again.");
  geomCmd->SetGuidance("This command must be applied if geometry has been modified");
  geomCmd->SetGuidance(" after the first initialization (or BeamOn).");
  geomCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  geomRebCmd = new G4UIcmdWithABool("/run/reinitializeGeometry",this);
  geomRebCmd->SetGuidance("Force geometry to be rebuilt once again.");
  geomRebCmd->SetGuidance("This command must be applied if the user needs his/her");
  geomRebCmd->SetGuidance(" detector construction to be reinvoked.");
  geomRebCmd->SetGuidance("/run/geometryModified is automatically issued with this command.");
  geomRebCmd->SetParameterName("destroyFirst",true);
  geomRebCmd->SetDefaultValue(false);
  geomRebCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  physCmd = new G4UIcmdWithoutParameter("/run/physicsModified",this);
  physCmd->SetGuidance("Force all physics tables recalculated again.");
  physCmd->SetGuidance("This command must be applied");
  physCmd->SetGuidance(" if physics process has been modified after the");
  physCmd->SetGuidance(" first initialization (or BeamOn).");
  physCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  constScoreCmd = new G4UIcmdWithoutParameter("/run/constructScoringWorlds",this);
  constScoreCmd->SetGuidance("Constrct scoring parallel world(s) if defined.");
  constScoreCmd->SetGuidance("This command is not mandatory, but automatically called when a run starts.");
  constScoreCmd->SetGuidance("But the user may use this to visualize the scoring world(s) before a run to start.");
  constScoreCmd->AvailableForStates(G4State_Idle);

  materialScanner = new G4MaterialScanner();
  
  randomDirectory = new G4UIdirectory("/random/");
  randomDirectory->SetGuidance("Random number status control commands.");

  seedCmd = new G4UIcmdWithAString("/random/setSeeds",this);
  seedCmd->SetGuidance("Initialize the random number generator with integer seed stream.");
  seedCmd->SetGuidance("Number of integers should be more than 1.");
  seedCmd->SetGuidance("Actual number of integers to be used depends on the individual random number engine.");
#ifdef G4MULTITHREADED
  seedCmd->SetGuidance("This command sets the seeds for the master thread.");
#endif
  seedCmd->SetParameterName("IntArray",false);
  seedCmd->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed);
  seedCmd->SetToBeBroadcasted(false);
  
  randDirCmd = new G4UIcmdWithAString("/random/setDirectoryName",this);
  randDirCmd->SetGuidance("Define the directory name of the rndm status files.");
  randDirCmd->SetGuidance("Directory will be created if it does not exist.");
  randDirCmd->SetParameterName("fileName",true);
  randDirCmd->SetDefaultValue("./");
  randDirCmd->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed);
  
  savingFlagCmd = new G4UIcmdWithABool("/random/setSavingFlag",this);
  savingFlagCmd->SetGuidance("The randomNumberStatus will be saved at :");
  savingFlagCmd->SetGuidance("begining of run (currentRun.rndm) and "
                             "begining of event (currentEvent.rndm) ");  
  savingFlagCmd->SetParameterName("flag",true);
  savingFlagCmd->SetDefaultValue(true);
  
  saveThisRunCmd = new G4UIcmdWithoutParameter("/random/saveThisRun",this);
  saveThisRunCmd->SetGuidance("copy currentRun.rndm to runXXX.rndm");
  saveThisRunCmd->AvailableForStates(G4State_Idle,G4State_GeomClosed,G4State_EventProc);
  
  saveThisEventCmd = new G4UIcmdWithoutParameter("/random/saveThisEvent",this);
  saveThisEventCmd->SetGuidance("copy currentEvent.rndm to runXXXevtYYY.rndm");
  saveThisEventCmd->AvailableForStates(G4State_EventProc);
          
  restoreRandCmd = new G4UIcmdWithAString("/random/resetEngineFrom",this);
  restoreRandCmd->SetGuidance("Reset the status of the rndm engine from a file.");
  restoreRandCmd->SetGuidance("See CLHEP manual for detail.");
  restoreRandCmd->SetGuidance("The engine status must be stored beforehand.");
  restoreRandCmd->SetGuidance("Directory of the status file should be set by"
                              " /random/setDirectoryName.");
  restoreRandCmd->SetParameterName("fileName",true);
  restoreRandCmd->SetDefaultValue("currentRun.rndm");
  restoreRandCmd->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed);
  restoreRandCmd->SetToBeBroadcasted(false);
    
  restoreRandCmdMT = new G4UIcmdWithABool("/random/resetEngineFromEachEvent",this);
  restoreRandCmdMT->SetGuidance("Reset the status of the rndm engine from a file at each event.");
  restoreRandCmdMT->SetGuidance("Note that the file must follow the following naming convention:");
  restoreRandCmdMT->SetGuidance("run{#1}evt{#2}.rndm ; where #1 is the run number and #2 is the event number.");
  restoreRandCmdMT->SetGuidance("For example to re-seed the first event of the first run the file should be called run0evt0.rndm.");
  restoreRandCmdMT->SetGuidance("If for a specific run/event the file is not found, the standard re-seeding strategy is used.");
  restoreRandCmdMT->SetGuidance("This command has meaning only in MT mode for strong reproducibility studies.");
  restoreRandCmdMT->SetGuidance("Directory of the status file should be set by"
                              	  " /random/setDirectoryName.");
  restoreRandCmdMT->SetDefaultValue(false);
  restoreRandCmdMT->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed);


  saveEachEventCmd = new G4UIcmdWithABool("/random/saveEachEventFlag",this);
  saveEachEventCmd->SetGuidance("Save random number status at beginning of each event.");
  saveEachEventCmd->SetGuidance("File name contains run and event numbers: runXXXevtYYY.rndm");
  saveEachEventCmd->SetParameterName("flag",true);
  saveEachEventCmd->SetDefaultValue(true);
  
  randEvtCmd = new G4UIcmdWithAnInteger("/run/storeRndmStatToEvent",this);
  randEvtCmd->SetGuidance("Flag to store rndm status to G4Event object.");
  randEvtCmd->SetGuidance(" flag = 0 : not store (default)");
  randEvtCmd->SetGuidance(" flag = 1 : status before primary particle generation is stored");
  randEvtCmd->SetGuidance(" flag = 2 : status before event processing (after primary particle generation) is stored");
  randEvtCmd->SetGuidance(" flag = 3 : both are stored");
  randEvtCmd->SetGuidance("Note: Some performance overhead may be seen by storing rndm status, in particular");
  randEvtCmd->SetGuidance("for the case of simplest geometry and small number of tracks per event.");
  randEvtCmd->SetParameterName("flag",true);
  randEvtCmd->SetDefaultValue(0);
  randEvtCmd->SetRange("flag>=0 && flag<3");
  randEvtCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
  procUICmds = new G4UIcmdWithoutParameter("/run/workersProcessCmds",this);
  procUICmds->SetToBeBroadcasted(false);
  procUICmds->SetGuidance("Force workers to process current stack of UI commands.");
  procUICmds->SetGuidance("This commands is meaningful only in MT mode.");
  procUICmds->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed);

}

G4RunMessenger::~G4RunMessenger()
{
  delete materialScanner;
  delete beamOnCmd;
  delete verboseCmd;
  delete printProgCmd;
  delete nThreadsCmd;
  delete maxThreadsCmd;
  delete pinAffinityCmd;
  delete evModCmd;
  delete optCmd;
  delete dumpRegCmd;
  delete dumpCoupleCmd;
  delete brkBoECmd;
  delete brkEoECmd;
  delete abortCmd;
  delete abortEventCmd;
  delete initCmd;
  delete geomCmd;
  delete geomRebCmd;
  delete physCmd;
  delete randEvtCmd;
  delete constScoreCmd;
  delete procUICmds;
  
  delete seedCmd;
  delete savingFlagCmd;
  delete saveThisRunCmd;
  delete saveThisEventCmd;
  delete restoreRandCmd;
  delete randomDirectory;
  delete saveEachEventCmd;
 
  delete randDirCmd;
  delete runDirectory;

  delete restoreRandCmdMT;
}

void G4RunMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==beamOnCmd )
  {
    G4int nev;
    G4int nst;
    const char* nv = (const char*)newValue;
    std::istringstream is(nv);
    is >> nev >> macroFileName >> nst;
    if(macroFileName=="***NULL***")
    { runManager->BeamOn(nev); }
    else
    { runManager->BeamOn(nev,macroFileName,nst); }
  }
  else if( command==verboseCmd )
  { runManager->SetVerboseLevel(verboseCmd->GetNewIntValue(newValue)); }
  else if( command == printProgCmd )
  { runManager->SetPrintProgress(printProgCmd->GetNewIntValue(newValue)); }
  else if( command==nThreadsCmd )
  {
    G4RunManager::RMType rmType = runManager->GetRunManagerType();
    if( rmType==G4RunManager::masterRM )
    {
      static_cast<G4MTRunManager*>(runManager)->SetNumberOfThreads(
       nThreadsCmd->GetNewIntValue(newValue));
    }
    else if ( rmType==G4RunManager::sequentialRM )
    {
      G4cout<<"*** /run/numberOfThreads command is issued in sequential mode."
            <<"\nCommand is ignored."<<G4endl;
    }
    else
    {
      G4Exception("G4RunMessenger::ApplyNewCommand","Run0901",FatalException,
      "/run/numberOfThreads command is issued to local thread.");
    }
  }
  else if( command==maxThreadsCmd)
  {
    G4RunManager::RMType rmType = runManager->GetRunManagerType();
    if( rmType==G4RunManager::masterRM )
    {
      static_cast<G4MTRunManager*>(runManager)->SetNumberOfThreads(
       G4Threading::G4GetNumberOfCores() );
    }
    else if ( rmType==G4RunManager::sequentialRM )
    {
      G4cout<<"*** /run/useMaximumLogicalCores command is issued in sequential mode."
            <<"\nCommand is ignored."<<G4endl;
    }
    else
    {
      G4Exception("G4RunMessenger::ApplyNewCommand","Run0901",FatalException,
      "/run/useMaximumLogicalCores command is issued to local thread.");
    }
  }
  else if ( command == pinAffinityCmd )
  {
    G4RunManager::RMType rmType = runManager->GetRunManagerType();
    if( rmType==G4RunManager::masterRM )
    {
      static_cast<G4MTRunManager*>(runManager)->SetPinAffinity(
    		  pinAffinityCmd->GetNewIntValue(newValue) );
    }
    else if ( rmType==G4RunManager::sequentialRM )
    {
      G4cout<<"*** /run/pinAffinity command is issued in sequential mode."
            <<"\nCommand is ignored."<<G4endl;
    }
    else
    {
      G4Exception("G4RunMessenger::ApplyNewCommand","Run0901",FatalException,
      "/run/pinAffinity command is issued to local thread.");
    }

  }
  else if( command==evModCmd)
  {
    G4RunManager::RMType rmType = runManager->GetRunManagerType();
    if( rmType==G4RunManager::masterRM )
    {
      G4int nevMod = 0;
      G4int sOnce = 0;
      const char* nv = (const char*)newValue;
      std::istringstream is(nv);
      is >> nevMod >> sOnce;
      static_cast<G4MTRunManager*>(runManager)->SetEventModulo(nevMod);
      G4MTRunManager::SetSeedOncePerCommunication(sOnce);
    }
    else if ( rmType==G4RunManager::sequentialRM )
    {
      G4cout<<"*** /run/eventModulo command is issued in sequential mode."
            <<"\nCommand is ignored."<<G4endl;
    }
    else
    {
      G4Exception("G4RunMessenger::ApplyNewCommand","Run0902",FatalException,
      "/run/eventModulo command is issued to local thread.");
    }
  }
  else if( command==dumpRegCmd )
  { 
    if(newValue=="**ALL**")
    { runManager->DumpRegion(); }
    else
    { runManager->DumpRegion(newValue); }
  }
  else if( command==dumpCoupleCmd)
  {
    G4ProductionCutsTable::GetProductionCutsTable()->DumpCouples();
  }
  else if( command==optCmd )
  { runManager->SetGeometryToBeOptimized(optCmd->GetNewBoolValue(newValue)); }
  else if( command==brkBoECmd )
  { G4UImanager::GetUIpointer()->SetPauseAtBeginOfEvent(brkBoECmd->GetNewBoolValue(newValue)); }
  else if( command==brkEoECmd )
  { G4UImanager::GetUIpointer()->SetPauseAtEndOfEvent(brkEoECmd->GetNewBoolValue(newValue)); }
  else if( command==abortCmd )
  { runManager->AbortRun(abortCmd->GetNewBoolValue(newValue)); }
  else if( command==abortEventCmd )
  { runManager->AbortEvent(); }
  else if( command==initCmd )
  { runManager->Initialize(); }
  else if( command==geomCmd )
  { runManager->GeometryHasBeenModified(false); }
  else if( command==geomRebCmd )
  { runManager->ReinitializeGeometry(geomRebCmd->GetNewBoolValue(newValue),false); }
  else if( command==physCmd )
  { runManager->PhysicsHasBeenModified(); }
  else if( command==seedCmd )
  {
    G4Tokenizer next(newValue);
    G4int idx=0;
    long seeds[100];
    G4String vl;
    while(!(vl=next()).isNull())
    { seeds[idx] = (long)(StoI(vl)); idx++; }
    if(idx<2)
    { G4cerr << "/random/setSeeds should have at least two integers. Command ignored." << G4endl; }
    else
    {
      seeds[idx] = 0;
      G4Random::setTheSeeds(seeds);
    }
  }
  else if( command==randDirCmd )
  { runManager->SetRandomNumberStoreDir(newValue); }
  else if( command==savingFlagCmd )
  { runManager->SetRandomNumberStore(savingFlagCmd->GetNewBoolValue(newValue)); }    
  else if( command==saveThisRunCmd )
  { runManager->rndmSaveThisRun(); }
  else if( command==saveThisEventCmd )
  { runManager->rndmSaveThisEvent(); }  
  else if( command==restoreRandCmd )
  { runManager->RestoreRandomNumberStatus(newValue); }
  else if( command==randEvtCmd )
  { runManager->StoreRandomNumberStatusToG4Event(randEvtCmd->GetNewIntValue(newValue)); }
  else if( command==saveEachEventCmd)
  { runManager->SetRandomNumberStorePerEvent(saveEachEventCmd->GetNewBoolValue(newValue)); }
  else if( command==constScoreCmd )
  { runManager->ConstructScoringWorlds(); }
  else if( command==restoreRandCmdMT)
  { runManager->RestoreRndmEachEvent(restoreRandCmdMT->GetNewBoolValue(newValue)); }
  else if ( command==procUICmds)
    {
      G4RunManager::RMType rmType = runManager->GetRunManagerType();
      if( rmType==G4RunManager::masterRM ) {
          auto rm = dynamic_cast<G4MTRunManager*>(runManager);
          if ( rm != nullptr ) {
              rm->RequestWorkersProcessCommandsStack();
          }
          else { G4Exception("G4RunManager::ApplyNewCommand","Run0128",FatalException,"/run/workersProcessCmds command issued on a non-G4MTRunManager class instance."); }
      } else if ( rmType==G4RunManager::sequentialRM ) {
          G4cout<<"*** /run/workersProcessCmds command is issued in sequential mode."
                <<"\nCommand is ignored."<<G4endl;
      } else {
          G4Exception("G4RunMessenger::ApplyNewCommand","Run0129",FatalException,
          "/run/workersProcessCmds command is issued to local thread.");
      }
    }

}

G4String G4RunMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  
  if( command==verboseCmd )
  { cv = verboseCmd->ConvertToString(runManager->GetVerboseLevel()); }
  else if( command == printProgCmd )
  { cv = printProgCmd->ConvertToString(runManager->GetPrintProgress()); }
  else if( command==randDirCmd )
  { cv = runManager->GetRandomNumberStoreDir(); }
  else if( command==randEvtCmd )
  { cv = randEvtCmd->ConvertToString(runManager->GetFlagRandomNumberStatusToG4Event()); }
  else if( command==nThreadsCmd )
  {
    G4RunManager::RMType rmType = runManager->GetRunManagerType();
    if( rmType==G4RunManager::masterRM )
    {
      cv = nThreadsCmd->ConvertToString(
       static_cast<G4MTRunManager*>(runManager)->GetNumberOfThreads());
    }
    else if ( rmType==G4RunManager::sequentialRM )
    {
      cv = "0";
    }
  }
  else if( command==evModCmd)
  {
    G4RunManager::RMType rmType = runManager->GetRunManagerType();
    if( rmType==G4RunManager::masterRM )
    {
      cv = evModCmd->ConvertToString(
       static_cast<G4MTRunManager*>(runManager)->GetEventModulo() )
       + " " + evModCmd->ConvertToString(G4MTRunManager::SeedOncePerCommunication());
    }
    else if ( rmType==G4RunManager::sequentialRM )
    { G4cout<<"*** /run/eventModulo command is valid only in MT mode."<<G4endl; }
  }
  
  return cv;
}


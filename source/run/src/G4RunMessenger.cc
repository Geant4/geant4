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
// $Id$
//

#include "G4RunMessenger.hh"
#include "G4RunManager.hh"
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

  verboseCmd = new G4UIcmdWithAnInteger("/run/verbose",this);
  verboseCmd->SetGuidance("Set the Verbose level of G4RunManager.");
  verboseCmd->SetGuidance(" 0 : Silent (default)");
  verboseCmd->SetGuidance(" 1 : Display main topics");
  verboseCmd->SetGuidance(" 2 : Display main topics and run summary");
  verboseCmd->SetParameterName("level",true);
  verboseCmd->SetDefaultValue(0);
  verboseCmd->SetRange("level >=0 && level <=2");

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
  geomCmd->SetGuidance("Force geometry to be closed again.");
  geomCmd->SetGuidance("This command must be applied");
  geomCmd->SetGuidance(" if geometry has been modified after the");
  geomCmd->SetGuidance(" first initialization (or BeamOn).");
  geomCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  physCmd = new G4UIcmdWithoutParameter("/run/physicsModified",this);
  physCmd->SetGuidance("Force all physics tables recalculated again.");
  physCmd->SetGuidance("This command must be applied");
  physCmd->SetGuidance(" if physics process has been modified after the");
  physCmd->SetGuidance(" first initialization (or BeamOn).");
  physCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  cutCmd = new G4UIcmdWithoutParameter("/run/cutoffModified",this);
  cutCmd->SetGuidance("/run/cutoffModified becomes obsolete.");
  cutCmd->SetGuidance("It is safe to remove invoking this command.");

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
  seedCmd->SetParameterName("IntArray",false);
  seedCmd->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed);
  
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

  //old commands for the rndm engine status handling
  //
  randDirOld = new G4UIcmdWithAString("/run/randomNumberStatusDirectory",this);
  randDirOld->SetGuidance("Define the directory name of the rndm status files.");
  randDirOld->SetGuidance("Directory must be creates before storing the files.");
  randDirOld->SetGuidance("OBSOLETE --- Please use commands in /random/ directory");
  randDirOld->SetParameterName("fileName",true);
  randDirOld->SetDefaultValue("./");
  randDirOld->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed);
  
  storeRandOld = new G4UIcmdWithAnInteger("/run/storeRandomNumberStatus",this);
  storeRandOld->SetGuidance("The randomNumberStatus will be saved at :");
  storeRandOld->SetGuidance("begining of run (currentRun.rndm) and "
                            "begining of event (currentEvent.rndm) ");  
  storeRandOld->SetGuidance("OBSOLETE --- Please use commands in /random/ directory");
  storeRandOld->SetParameterName("flag",true);
  storeRandOld->SetDefaultValue(1);
          
  restoreRandOld = new G4UIcmdWithAString("/run/restoreRandomNumberStatus",this);
  restoreRandOld->SetGuidance("Reset the status of the rndm engine from a file.");
  restoreRandOld->SetGuidance("See CLHEP manual for detail.");
  restoreRandOld->SetGuidance("The engine status must be stored beforehand.");
  restoreRandOld->SetGuidance("Directory of the status file should be set by"
                              " /random/setDirectoryName.");
  restoreRandOld->SetGuidance("OBSOLETE --- Please use commands in /random/ directory");
  restoreRandOld->SetParameterName("fileName",true);
  restoreRandOld->SetDefaultValue("currentRun.rndm");
  restoreRandOld->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed);  
}

G4RunMessenger::~G4RunMessenger()
{
  delete materialScanner;
  delete beamOnCmd;
  delete verboseCmd;
  delete optCmd;
  delete dumpRegCmd;
  delete dumpCoupleCmd;
  delete brkBoECmd;
  delete brkEoECmd;
  delete abortCmd;
  delete abortEventCmd;
  delete initCmd;
  delete geomCmd;
  delete physCmd;
  delete cutCmd;
  delete randEvtCmd;
  delete randDirOld; delete storeRandOld; delete restoreRandOld; 
  delete constScoreCmd;
  delete runDirectory;
  
  delete randDirCmd;
  delete seedCmd;
  delete savingFlagCmd;
  delete saveThisRunCmd;
  delete saveThisEventCmd;
  delete restoreRandCmd;
  delete randomDirectory;
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
  { runManager->GeometryHasBeenModified(); }
  else if( command==physCmd )
  { runManager->PhysicsHasBeenModified(); }
  else if( command==cutCmd )
  { runManager->CutOffHasBeenModified(); }
 
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
      CLHEP::HepRandom::setTheSeeds(seeds);
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
  
  else if( command==randDirOld )
  {G4cout << "warning: deprecated command. Use /random/setDirectoryName"
          << G4endl; 
  // runManager->SetRandomNumberStoreDir(newValue);
  }
  else if( command==storeRandOld )
  {G4cout << "warning: deprecated command. Use /random/setSavingFlag"
          << G4endl;
   // G4int frequency = storeRandOld->GetNewIntValue(newValue);
   // G4bool flag = false;
   // if(frequency != 0) flag = true;	     
   // runManager->SetRandomNumberStore(flag);
  }    
  else if( command==restoreRandOld )
  {G4cout << "warning: deprecated command. Use /random/resetEngineFrom"
           << G4endl;  
   // runManager->RestoreRandomNumberStatus(newValue);
  }  
  else if( command==constScoreCmd )
  { runManager->ConstructScoringWorlds(); }

}

G4String G4RunMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  
  if( command==verboseCmd )
  { cv = verboseCmd->ConvertToString(runManager->GetVerboseLevel()); }
  else if( command==randDirCmd )
  { cv = runManager->GetRandomNumberStoreDir(); }
  else if( command==randEvtCmd )
  { cv = randEvtCmd->ConvertToString(runManager->GetFlagRandomNumberStatusToG4Event()); }
  
  return cv;
}


//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4RunMessenger.cc,v 1.10 2001-11-23 16:20:30 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4ios.hh"
#include "g4std/strstream"

G4RunMessenger::G4RunMessenger(G4RunManager * runMgr)
:runManager(runMgr)
{
  runDirectory = new G4UIdirectory("/run/");
  runDirectory->SetGuidance("Run control commands.");

  initCmd = new G4UIcmdWithoutParameter("/run/initialize",this);
  initCmd->SetGuidance("Initialize G4 kernel.");
  initCmd->AvailableForStates(PreInit,Idle),

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
  beamOnCmd->AvailableForStates(PreInit,Idle);
  G4UIparameter* p1 = new G4UIparameter("numberOfEvent",'i',true);
  p1->SetDefaultValue(1);
  p1->SetParameterRange("numberOfEvent > 0");
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

  optCmd = new G4UIcmdWithABool("/run/optimizeGeometry",this);
  optCmd->SetGuidance("Set the optimization flag for geometry.");
  optCmd->SetGuidance("If it is set to TRUE, G4GeometryManager will optimize");
  optCmd->SetGuidance("the geometry definitions.");
  optCmd->SetGuidance("GEANT4 is initialized with this flag as TRUE.");
  optCmd->SetParameterName("optimizeFlag",true);
  optCmd->SetDefaultValue(true);
  optCmd->AvailableForStates(PreInit,Idle);

  brkBoECmd = new G4UIcmdWithABool("/run/breakAtBeginOfEvent",this);
  brkBoECmd->SetGuidance("Set a break point at the begining of every event.");
  brkBoECmd->SetParameterName("flag",true);
  brkBoECmd->SetDefaultValue(true);
  
  brkEoECmd = new G4UIcmdWithABool("/run/breakAtEndOfEvent",this);
  brkEoECmd->SetGuidance("Set a break point at the end of every event.");
  brkEoECmd->SetParameterName("flag",true);
  brkEoECmd->SetDefaultValue(true);
  
  abortCmd = new G4UIcmdWithoutParameter("/run/abort",this);
  abortCmd->SetGuidance("Abort current run processing.");
  abortCmd->AvailableForStates(GeomClosed,EventProc);

  geomCmd = new G4UIcmdWithoutParameter("/run/geometryModified",this);
  geomCmd->SetGuidance("Force geometry to be closed again.");
  geomCmd->SetGuidance("This command must be applied");
  geomCmd->SetGuidance(" if geometry has been modified after the");
  geomCmd->SetGuidance(" first initialization (or BeamOn).");
  geomCmd->AvailableForStates(PreInit,Idle),

  cutCmd = new G4UIcmdWithoutParameter("/run/cutoffModified",this);
  cutCmd->SetGuidance("Force closssection tables to be calculated again.");
  cutCmd->SetGuidance("This command must be applied");
  cutCmd->SetGuidance(" if cutoff value(s) have been modified after the");
  cutCmd->SetGuidance(" first initialization (or BeamOn).");
  cutCmd->AvailableForStates(PreInit,Idle);
  
  randomDirectory = new G4UIdirectory("/random/");
  randomDirectory->SetGuidance("Random number status control commands.");
  
  randDirCmd = new G4UIcmdWithAString("/random/setDirectoryName",this);
  randDirCmd->SetGuidance("Define the directory name of the rndm status files.");
  randDirCmd->SetGuidance("Directory must be creates before storing the files.");
  randDirCmd->SetParameterName("fileName",true);
  randDirCmd->SetDefaultValue("./");
  randDirCmd->AvailableForStates(PreInit,Idle,GeomClosed);
  
  savingFlagCmd = new G4UIcmdWithABool("/random/setSavingFlag",this);
  savingFlagCmd->SetGuidance("The randomNumberStatus will be saved at :");
  savingFlagCmd->SetGuidance("begining of run (currentRun.rndm) and "
                             "begining of event (currentEvent.rndm) ");  
  savingFlagCmd->SetParameterName("flag",true);
  savingFlagCmd->SetDefaultValue(true);
  
  saveThisRunCmd = new G4UIcmdWithoutParameter("/random/saveThisRun",this);
  saveThisRunCmd->SetGuidance("copy currentRun.rndm to runXXX.rndm");
  saveThisRunCmd->AvailableForStates(Idle,GeomClosed,EventProc);
  
  saveThisEventCmd = new G4UIcmdWithoutParameter("/random/saveThisEvent",this);
  saveThisEventCmd->SetGuidance("copy currentEvent.rndm to runXXXevtYYY.rndm");
  saveThisEventCmd->AvailableForStates(EventProc);
          
  restoreRandCmd = new G4UIcmdWithAString("/random/resetEngineFrom",this);
  restoreRandCmd->SetGuidance("Reset the status of the rndm engine from a file.");
  restoreRandCmd->SetGuidance("See CLHEP manual for detail.");
  restoreRandCmd->SetGuidance("The engine status must be stored beforehand.");
  restoreRandCmd->SetGuidance("Directory of the status file should be set by"
                              " /random/setDirectoryName.");
  restoreRandCmd->SetParameterName("fileName",true);
  restoreRandCmd->SetDefaultValue("currentRun.rndm");
  restoreRandCmd->AvailableForStates(PreInit,Idle,GeomClosed);
  
  //old commands for the rndm engine status handling
  //
  randDirOld = new G4UIcmdWithAString("/run/randomNumberStatusDirectory",this);
  randDirOld->SetGuidance("Define the directory name of the rndm status files.");
  randDirOld->SetGuidance("Directory must be creates before storing the files.");
  randDirOld->SetParameterName("fileName",true);
  randDirOld->SetDefaultValue("./");
  randDirOld->AvailableForStates(PreInit,Idle,GeomClosed);
  
  storeRandOld = new G4UIcmdWithAnInteger("/run/storeRandomNumberStatus",this);
  storeRandOld->SetGuidance("The randomNumberStatus will be saved at :");
  storeRandOld->SetGuidance("begining of run (currentRun.rndm) and "
                            "begining of event (currentEvent.rndm) ");  
  storeRandOld->SetParameterName("flag",true);
  storeRandOld->SetDefaultValue(1);
          
  restoreRandOld = new G4UIcmdWithAString("/run/restoreRandomNumberStatus",this);
  restoreRandOld->SetGuidance("Reset the status of the rndm engine from a file.");
  restoreRandOld->SetGuidance("See CLHEP manual for detail.");
  restoreRandOld->SetGuidance("The engine status must be stored beforehand.");
  restoreRandOld->SetGuidance("Directory of the status file should be set by"
                              " /random/setDirectoryName.");
  restoreRandOld->SetParameterName("fileName",true);
  restoreRandOld->SetDefaultValue("currentRun.rndm");
  restoreRandOld->AvailableForStates(PreInit,Idle,GeomClosed);  

}

G4RunMessenger::~G4RunMessenger()
{
  delete beamOnCmd;
  delete verboseCmd;
  delete optCmd;
  delete brkBoECmd;
  delete brkEoECmd;
  delete abortCmd;
  delete initCmd;
  delete geomCmd;
  delete cutCmd;
  delete randDirOld; delete storeRandOld; delete restoreRandOld; 
  delete runDirectory;
  
  delete randDirCmd;
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
    G4int ns;
    const char* nv = (const char*)newValue;
    G4std::istrstream is((char*)nv);
    is >> nev >> macroFileName >> ns;
    if(macroFileName=="***NULL***")
    { runManager->BeamOn(nev); }
    else
    { runManager->BeamOn(nev,macroFileName,ns); }
  }
  else if( command==verboseCmd )
  { runManager->SetVerboseLevel(verboseCmd->GetNewIntValue(newValue)); }
  else if( command==optCmd )
  { runManager->SetGeometryToBeOptimized(optCmd->GetNewBoolValue(newValue)); }
  else if( command==brkBoECmd )
  { G4UImanager::GetUIpointer()->SetPauseAtBeginOfEvent(brkBoECmd->GetNewBoolValue(newValue)); }
  else if( command==brkEoECmd )
  { G4UImanager::GetUIpointer()->SetPauseAtEndOfEvent(brkEoECmd->GetNewBoolValue(newValue)); }
  else if( command==abortCmd )
  { runManager->AbortRun(); }
  else if( command==initCmd )
  { runManager->Initialize(); }
  else if( command==geomCmd )
  { runManager->GeometryHasBeenModified(); }
  else if( command==cutCmd )
  { runManager->CutOffHasBeenModified(); }
  
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
  
  else if( command==randDirOld )
  {G4cout << "warning: deprecated command. Use /random/setDirectoryName"
          << G4endl; 
   runManager->SetRandomNumberStoreDir(newValue); }
  else if( command==storeRandOld )
  {G4cout << "warning: deprecated command. Use /random/setSavingFlag"
          << G4endl;
   G4int frequency = storeRandOld->GetNewIntValue(newValue);
   G4bool flag = false;
   if(frequency != 0) flag = true;	     
   runManager->SetRandomNumberStore(flag); }    
  else if( command==restoreRandOld )
  {G4cout << "warning: deprecated command. Use /random/resetEngineFrom"
           << G4endl;  
   runManager->RestoreRandomNumberStatus(newValue); }  
}

G4String G4RunMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  
  if( command==verboseCmd )
  { cv = verboseCmd->ConvertToString(runManager->GetVerboseLevel()); }
  else if( command==randDirCmd )
  { cv = runManager->GetRandomNumberStoreDir(); }
  
  return cv;
}


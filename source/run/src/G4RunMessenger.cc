// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RunMessenger.cc,v 1.4 1999-07-25 15:58:47 asaim Exp $
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
#include "G4ios.hh"
#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif

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

  storeRandCmd = new G4UIcmdWithAnInteger("/run/storeRandomNumberStatus",this);
  storeRandCmd->SetGuidance("Store the status of the random number engine to a file.");
  storeRandCmd->SetGuidance("See CLHEP manual for detail.");
  storeRandCmd->SetGuidance("Frequency : 0 - not stored");
  storeRandCmd->SetGuidance("            1 - begining of each run");
  storeRandCmd->SetGuidance("           -1 - begining of each run (file overwritten)");
  storeRandCmd->SetGuidance("            2 - begining of each event before generating primaries");
  storeRandCmd->SetGuidance("           -2 - begining of each event before generating primaries (file overwritten)");
  storeRandCmd->SetGuidance("Stored status can be restored by ");
  storeRandCmd->SetGuidance("  /run/restoreRandomNumberStatus command.");
  storeRandCmd->SetGuidance("In case Frequency is negative, a output file (file name RandEngine.stat)");
  storeRandCmd->SetGuidance("  is overwitten every time.");
  storeRandCmd->SetGuidance("In case Frequency is 1, file names are RandEngineRxxx.stat,");
  storeRandCmd->SetGuidance("  where xxx is the run number.");
  storeRandCmd->SetGuidance("In case Frequency is 2, file names are RandEngineRxxxEyyy.stat,");
  storeRandCmd->SetGuidance("  where xxx is the run number and yyy is the event number.");
  storeRandCmd->AvailableForStates(PreInit,Idle);
  storeRandCmd->SetParameterName("frequency",true);
  storeRandCmd->SetDefaultValue(1);
  storeRandCmd->SetRange("frequency >=-2 && frequency <=2");

  restoreRandCmd = new G4UIcmdWithAString("/run/restoreRandomNumberStatus",this);
  restoreRandCmd->SetGuidance("Restore the status of the random number engine from a file.");
  restoreRandCmd->SetGuidance("See CLHEP manual for detail.");
  restoreRandCmd->SetGuidance("The engine status must be stored beforehand.");
  restoreRandCmd->SetParameterName("fileName",true);
  restoreRandCmd->SetDefaultValue("RandEngine.stat");
  restoreRandCmd->AvailableForStates(PreInit,Idle,GeomClosed);

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
  delete storeRandCmd;
  delete restoreRandCmd;
  delete runDirectory;
}

void G4RunMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==beamOnCmd )
  {
    G4int nev;
    G4int ns;
    const char* nv = (const char*)newValue;
    istrstream is((char*)nv);
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
  { runManager->SetPauseAtBeginOfEvent(brkBoECmd->GetNewBoolValue(newValue)); }
  else if( command==brkEoECmd )
  { runManager->SetPauseAtEndOfEvent(brkEoECmd->GetNewBoolValue(newValue)); }
  else if( command==abortCmd )
  { runManager->AbortRun(); }
  else if( command==initCmd )
  { runManager->Initialize(); }
  else if( command==geomCmd )
  { runManager->GeometryHasBeenModified(); }
  else if( command==cutCmd )
  { runManager->CutOffHasBeenModified(); }
  else if( command==storeRandCmd )
  { runManager->SetRandomNumberStore(storeRandCmd->GetNewIntValue(newValue)); }
  else if( command==restoreRandCmd )
  { runManager->RestoreRandomNumberStatus(newValue); }
}

G4String G4RunMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  
  if( command==verboseCmd )
  { cv = verboseCmd->ConvertToString(runManager->GetVerboseLevel()); }
  else if( command==storeRandCmd )
  { cv = storeRandCmd->ConvertToString(runManager->GetRandomNumberStore()); }
  
  return cv;
}


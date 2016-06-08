// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PersistencyMessenger.cc,v 1.7 1999/11/26 16:50:21 morita Exp $
// GEANT4 tag $Name: geant4-02-00 $
//

#include "G4PersistencyMessenger.hh"
#include "G4PersistencyManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4ios.hh"

G4PersistencyMessenger::G4PersistencyMessenger(G4PersistencyManager* persistencyMgr)
: persistencyManager(persistencyMgr)
{
  persistencyDirectory = new G4UIdirectory("/db/");
  persistencyDirectory->SetGuidance("Database control commands.");

  verboseCmd = new G4UIcmdWithAnInteger("/db/verbose",this);
  verboseCmd->SetGuidance("Set the verbose level of G4PersistencyManager.");
  verboseCmd->SetGuidance(" 0 : Silent (default)");
  verboseCmd->SetGuidance(" 1 : Display main topics");
  verboseCmd->SetGuidance(" 2 : Display event-level topics");
  verboseCmd->SetGuidance(" 3 : Display debug information");
  verboseCmd->SetParameterName("level",true);
  verboseCmd->SetDefaultValue(0);
  verboseCmd->SetRange("level >=0 && level <=3");

  runDbCmd = new G4UIcmdWithAString("/db/run",this);
  runDbCmd->SetGuidance("Set the name of Run Database.");
  runDbCmd->SetParameterName("fileName",true);
  runDbCmd->SetDefaultValue("Runs");
  runDbCmd->AvailableForStates(PreInit,Idle);

  eventDbCmd = new G4UIcmdWithAString("/db/event",this);
  eventDbCmd->SetGuidance("Set the name of Event Database.");
  eventDbCmd->SetParameterName("fileName",true);
  eventDbCmd->SetDefaultValue("Events");
  eventDbCmd->AvailableForStates(PreInit,Idle);

  geomDbCmd = new G4UIcmdWithAString("/db/geometry",this);
  geomDbCmd->SetGuidance("Set the name of Geometry Database.");
  geomDbCmd->SetParameterName("fileName",true);
  geomDbCmd->SetDefaultValue("Geometry");
  geomDbCmd->AvailableForStates(PreInit,Idle);
}

G4PersistencyMessenger::~G4PersistencyMessenger()
{
  delete verboseCmd;
  delete runDbCmd;
  delete eventDbCmd;
  delete geomDbCmd;
  delete persistencyDirectory;
}

void G4PersistencyMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command==verboseCmd )
  { persistencyManager->SetVerboseLevel(verboseCmd->GetNewIntValue(newValue)); }
  else if( command==runDbCmd )
  { persistencyManager->SelectDB(kRunDB, newValue, true); }
  else if( command==eventDbCmd )
  { persistencyManager->SelectDB(kEventDB, newValue, true); }
  else if( command==geomDbCmd )
  { persistencyManager->SelectDB(kGeomDB, newValue, true); }
}

G4String G4PersistencyMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;

  if( command==verboseCmd )
  { cv = verboseCmd->ConvertToString(persistencyManager->GetVerboseLevel()); }
  else if( command==runDbCmd )
  { cv = persistencyManager->DBName(kRunDB); }
  else if( command==eventDbCmd )
  { cv = persistencyManager->DBName(kEventDB); }
  else if( command==geomDbCmd )
  { cv = persistencyManager->DBName(kGeomDB); }

  return cv;
}


// Rich advanced example for Geant4
// RichTbPrimaryGeneratorMessenger.cc for Rich of LHCb
// History:
// Created: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include "RichTbPrimaryGeneratorMessenger.hh"

#include "RichTbPrimaryGeneratorAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RichTbPrimaryGeneratorMessenger::RichTbPrimaryGeneratorMessenger(RichTbPrimaryGeneratorAction* RichTbGun)
:RichTbAction(RichTbGun)
{ 
  RndmCmd = new G4UIcmdWithAString("/gun/random",this);
  RndmCmd->SetGuidance("Shoot randomly the incident particle.");
  RndmCmd->SetGuidance("  Choice : on, off(default)");
  RndmCmd->SetParameterName("choice",true);
  RndmCmd->SetDefaultValue("off");
  RndmCmd->SetCandidates("on off");
  RndmCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  setxvertexCmd = new G4UIcmdWithADoubleAndUnit("/gun/xvertex",this);
  setxvertexCmd->SetGuidance(" Set x coord. of the primary vertex.");
  setxvertexCmd->SetParameterName("xv",true);
  setxvertexCmd->SetDefaultValue(0.0*mm) ; 
  
  setyvertexCmd = new G4UIcmdWithADoubleAndUnit("/gun/yvertex",this);
  setyvertexCmd->SetGuidance(" Set y coord. of the primary vertex.");
  setyvertexCmd->SetParameterName("yv",true);
  setyvertexCmd->SetDefaultValue(0.0*mm) ; 
  
  setzvertexCmd = new G4UIcmdWithADoubleAndUnit("/gun/zvertex",this);
  setzvertexCmd->SetGuidance(" Set z coord. of the primary vertex.");
  setzvertexCmd->SetParameterName("zv",true);
  setzvertexCmd->SetDefaultValue(0.0*mm) ; 
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RichTbPrimaryGeneratorMessenger::~RichTbPrimaryGeneratorMessenger()
{
  delete RndmCmd;
  delete setxvertexCmd;
  delete setyvertexCmd;
  delete setzvertexCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RichTbPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if( command == RndmCmd )
   { RichTbAction->SetRndmFlag(newValue);}
  if( command == setxvertexCmd)
   { RichTbAction->Setxvertex(setxvertexCmd->GetNewDoubleValue(newValue));}
  if( command == setyvertexCmd)
   { RichTbAction->Setyvertex(setyvertexCmd->GetNewDoubleValue(newValue));}
  if( command == setzvertexCmd)
   { RichTbAction->Setzvertex(setzvertexCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


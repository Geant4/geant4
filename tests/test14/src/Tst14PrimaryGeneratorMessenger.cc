// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst14PrimaryGeneratorMessenger.cc,v 1.1 1999-05-29 14:12:12 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Tst14PrimaryGeneratorMessenger.hh"

#include "Tst14PrimaryGeneratorAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst14PrimaryGeneratorMessenger::Tst14PrimaryGeneratorMessenger(Tst14PrimaryGeneratorAction* Tst14Gun)
:Tst14Action(Tst14Gun)
{ 
  RndmCmd = new G4UIcmdWithAString("/gun/random",this);
  RndmCmd->SetGuidance("Shoot randomly the incident particle.");
  RndmCmd->SetGuidance("  Choice : on, off(default)");
  RndmCmd->SetParameterName("choice",true);
  RndmCmd->SetDefaultValue("off");
  RndmCmd->SetCandidates("on off");
  RndmCmd->AvailableForStates(PreInit,Idle);
 
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

Tst14PrimaryGeneratorMessenger::~Tst14PrimaryGeneratorMessenger()
{
  delete RndmCmd;
  delete setxvertexCmd;
  delete setyvertexCmd;
  delete setzvertexCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst14PrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if( command == RndmCmd )
   { Tst14Action->SetRndmFlag(newValue);}
  if( command == setxvertexCmd)
   { Tst14Action->Setxvertex(setxvertexCmd->GetNewDoubleValue(newValue));}
  if( command == setyvertexCmd)
   { Tst14Action->Setyvertex(setyvertexCmd->GetNewDoubleValue(newValue));}
  if( command == setzvertexCmd)
   { Tst14Action->Setzvertex(setzvertexCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


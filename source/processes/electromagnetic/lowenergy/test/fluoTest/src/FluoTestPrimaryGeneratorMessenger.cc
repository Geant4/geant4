// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: FluoTestPrimaryGeneratorMessenger.cc,v 1.4 2001-10-25 16:35:47 elena Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "FluoTestPrimaryGeneratorMessenger.hh"

#include "FluoTestPrimaryGeneratorAction.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestPrimaryGeneratorMessenger::FluoTestPrimaryGeneratorMessenger(FluoTestPrimaryGeneratorAction* FluoTestGun)
:FluoTestAction(FluoTestGun)
{ 
  RndmCmd = new G4UIcmdWithAString("/gun/random",this);
  RndmCmd->SetGuidance("Shoot randomly the incident particle.");
  RndmCmd->SetGuidance("  Choice : on(default), off");
  RndmCmd->SetParameterName("choice",true);
  RndmCmd->SetDefaultValue("on");
  RndmCmd->SetCandidates("on off");
  RndmCmd->AvailableForStates(PreInit,Idle);

 RndmPart = new G4UIcmdWithAString("/gun/randomPart",this);
  RndmPart->SetGuidance("Shoot randomly the incident particle.");
  RndmPart->SetGuidance("  Choice : on(default), off");
  RndmPart->SetParameterName("choice",true);
  RndmPart->SetDefaultValue("on");
  RndmPart->SetCandidates("on off");
  RndmPart->AvailableForStates(PreInit,Idle);

  RndmVert = new G4UIcmdWithAString("/gun/randomVert",this);
  RndmVert->SetGuidance("Shoot randomly the incident particle.");
  RndmVert->SetGuidance("  Choice : on(default), off");
  RndmVert->SetParameterName("choice",true);
  RndmVert->SetDefaultValue("on");
  RndmVert->SetCandidates("on off");
  RndmVert->AvailableForStates(PreInit,Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestPrimaryGeneratorMessenger::~FluoTestPrimaryGeneratorMessenger()
{
  delete RndmCmd;
  delete RndmPart;
  delete  RndmVert;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FluoTestPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if( command == RndmCmd )
   { FluoTestAction->SetRndmFlag(newValue);}
 if( command == RndmPart )
   { FluoTestAction->SetRndmPart(newValue);}
 if( command == RndmVert )
   { FluoTestAction->SetRndmVert(newValue);}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em3PrimaryGeneratorMessenger.cc,v 1.3 2001-04-13 13:17:32 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em3PrimaryGeneratorMessenger.hh"

#include "Em3PrimaryGeneratorAction.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em3PrimaryGeneratorMessenger::Em3PrimaryGeneratorMessenger(Em3PrimaryGeneratorAction* Em3Gun)
:Em3Action(Em3Gun)
{ 
  DefaultCmd = new G4UIcmdWithoutParameter("/gun/setDefault",this);
  DefaultCmd->SetGuidance("set/reset the kinematic defined in PrimaryGenerator");
  DefaultCmd->AvailableForStates(PreInit,Idle);
  
  RndmCmd = new G4UIcmdWithADouble("/gun/rndm",this);
  RndmCmd->SetGuidance("random lateral extension on the beam");
  RndmCmd->SetGuidance("in fraction of 0.5*sizeYZ");
  RndmCmd->SetParameterName("rBeam",false);
  RndmCmd->SetRange("rBeam>=0.&&rBeam<=1.");
  RndmCmd->AvailableForStates(Idle);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em3PrimaryGeneratorMessenger::~Em3PrimaryGeneratorMessenger()
{
  delete DefaultCmd;
  delete RndmCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em3PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == DefaultCmd )
   { Em3Action->SetDefaultKinematic();}
   
  if( command == RndmCmd )
   { Em3Action->SetRndmBeam(RndmCmd->GetNewDoubleValue(newValue));}   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: XrayFluoPrimaryGeneratorMessenger.cc,v 1.1 2001-11-27 14:59:35 elena Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "XrayFluoPrimaryGeneratorMessenger.hh"

#include "XrayFluoPrimaryGeneratorAction.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoPrimaryGeneratorMessenger::XrayFluoPrimaryGeneratorMessenger(XrayFluoPrimaryGeneratorAction* XrayFluoGun)
:XrayFluoAction(XrayFluoGun)
{ 
  RndmCmd = new G4UIcmdWithAString("/gun/random",this);
  RndmCmd->SetGuidance("Shoot randomly the incident particle.");
  RndmCmd->SetGuidance("  Choice : on(default), off");
  RndmCmd->SetParameterName("choice",true);
  RndmCmd->SetDefaultValue("on");
  RndmCmd->SetCandidates("on off");
  RndmCmd->AvailableForStates(PreInit,Idle);

  RndmVert = new G4UIcmdWithAString("/gun/randomVert",this);
  RndmVert->SetGuidance("Shoot randomly the incident particle.");
  RndmVert->SetGuidance("  Choice : on(default), off");
  RndmVert->SetParameterName("choice",true);
  RndmVert->SetDefaultValue("on");
  RndmVert->SetCandidates("on off");
  RndmVert->AvailableForStates(PreInit,Idle);

  spectrum = new G4UIcmdWithAString("/gun/spectrum",this);
  spectrum->SetGuidance("Shoot the incident particle with a certain energy spectrum.");
  spectrum->SetGuidance("  Choice : on(default), off");
  spectrum->SetParameterName("choice",true);
  spectrum->SetDefaultValue("on");
  spectrum->SetCandidates("on off");
  spectrum->AvailableForStates(PreInit,Idle);

  isoVert = new G4UIcmdWithAString("/gun/isoVert",this);
  isoVert->SetGuidance("Shoot the incident particle from an isotrofic direction.");
  isoVert->SetGuidance("  Choice : on(default), off");
  isoVert->SetParameterName("choice",true);
  isoVert->SetDefaultValue("on");
  isoVert->SetCandidates("on off");
  isoVert->AvailableForStates(PreInit,Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoPrimaryGeneratorMessenger::~XrayFluoPrimaryGeneratorMessenger()
{
  delete RndmCmd;
  delete  RndmVert;
  delete spectrum;
  delete isoVert;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if( command == RndmCmd )
   { XrayFluoAction->SetRndmFlag(newValue);}
 if( command == RndmVert )
   { XrayFluoAction->SetRndmVert(newValue);}
 if( command == spectrum )
   { XrayFluoAction->SetSpectrum(newValue);} 
 if( command == isoVert )
   { XrayFluoAction->SetIsoVert(newValue);}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


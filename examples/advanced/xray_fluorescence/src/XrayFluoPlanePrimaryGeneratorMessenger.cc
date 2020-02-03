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
//
// Author: Alfonso Mantero (Alfonso.mantero@ge.infn.it)
//
// History:
// -----------
// 02 Sep 2003  Alfonso Mantero created
//
// -------------------------------------------------------------------

#include "XrayFluoPlanePrimaryGeneratorMessenger.hh"
#include "XrayFluoPlanePrimaryGeneratorAction.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoPlanePrimaryGeneratorMessenger::XrayFluoPlanePrimaryGeneratorMessenger(XrayFluoPlanePrimaryGeneratorAction* XrayFluoGun)
:XrayFluoAction(XrayFluoGun)
{ 
  RndmCmd = new G4UIcmdWithAString("/gun/random",this);
  RndmCmd->SetGuidance("Shoot particles from a point-like source at infinity, witth an incidence angle with the plane of 45 deg");
  RndmCmd->SetGuidance("  Choice : on(default), off");
  RndmCmd->SetParameterName("choice",true);
  RndmCmd->SetDefaultValue("on");
  RndmCmd->SetCandidates("on off");
  RndmCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  RndmVert = new G4UIcmdWithAString("/gun/beam",this);
  RndmVert->SetGuidance("Creates a round beam of particles of 0.5mm in diameter.");
  RndmVert->SetGuidance("  Choice : on(default), off");
  RndmVert->SetParameterName("choice",true);
  RndmVert->SetDefaultValue("on");
  RndmVert->SetCandidates("on off");
  RndmVert->AvailableForStates(G4State_PreInit,G4State_Idle);

  spectrum = new G4UIcmdWithAString("/gun/spectrum",this);
  spectrum->SetGuidance("Shoot the incident particle with a certain energy spectrum.");
  spectrum->SetGuidance("  Choice : on(default), off");
  spectrum->SetParameterName("choice",true);
  spectrum->SetDefaultValue("on");
  spectrum->SetCandidates("on off");
  spectrum->AvailableForStates(G4State_PreInit,G4State_Idle);

  isoVert = new G4UIcmdWithAString("/gun/isoVert",this);
  isoVert->SetGuidance("Shoot the incident particle from an isotropic direction.");
  isoVert->SetGuidance("  Choice : on(default), off");
  isoVert->SetParameterName("choice",true);
  isoVert->SetDefaultValue("on");
  isoVert->SetCandidates("on off");
  isoVert->AvailableForStates(G4State_PreInit,G4State_Idle);

  G4cout << "XrayFluoPlanePrimaryGeneratorMessenger created" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoPlanePrimaryGeneratorMessenger::~XrayFluoPlanePrimaryGeneratorMessenger()
{
  delete RndmCmd;
  delete  RndmVert;
  delete spectrum;
  delete isoVert;

  G4cout << "XrayFluoPlanePrimaryGeneratorMessenger deleted" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPlanePrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if( command == RndmCmd )
    { 
      XrayFluoAction->SetRndmFlag(newValue);
      XrayFluoAction->SetRndmVert("off");
      XrayFluoAction->SetIsoVert("off");
    }
  
  if( command == RndmVert )
    { 
      XrayFluoAction->SetRndmVert(newValue);
      XrayFluoAction->SetIsoVert("off");
      XrayFluoAction->SetRndmFlag("off");
    }
  if( command == spectrum )
    { XrayFluoAction->SetSpectrum(newValue);} 
  
  if( command == isoVert )
    { 
      XrayFluoAction->SetIsoVert(newValue);
      XrayFluoAction->SetRndmFlag("off");
      XrayFluoAction->SetRndmVert("off");
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


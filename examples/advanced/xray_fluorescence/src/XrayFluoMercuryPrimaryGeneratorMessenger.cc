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
// 19 Sep 2003  Alfonso Mantero created
//
// -------------------------------------------------------------------

#include "XrayFluoMercuryPrimaryGeneratorMessenger.hh"
#include "XrayFluoMercuryPrimaryGeneratorAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoMercuryPrimaryGeneratorMessenger::XrayFluoMercuryPrimaryGeneratorMessenger(XrayFluoMercuryPrimaryGeneratorAction* XrayFluoGun)
:XrayFluoAction(XrayFluoGun)
{ 
  spectrum = new G4UIcmdWithAString("/gun/spectrum",this);
  spectrum->SetGuidance("Shoot the incident particle with a certain energy spectrum.");
  spectrum->SetGuidance("  Choice : on(default), off");
  spectrum->SetParameterName("choice",true);
  spectrum->SetDefaultValue("on");
  spectrum->SetCandidates("on off");
  spectrum->AvailableForStates(G4State_PreInit,G4State_Idle);

  globalCmd = new G4UIcmdWithABool("/gun/globalIllumunation",this);
  globalCmd->SetGuidance("Illuminate the entire Mercury globe, not only the area seen by the optics");
  globalCmd->SetGuidance("Choice : true, 1, false, 0");
  globalCmd->SetParameterName("illuminaton",true);
  globalCmd->SetDefaultValue("false");
  globalCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  G4cout << "XrayFluoMercuryPrimaryGeneratorMessenger created" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoMercuryPrimaryGeneratorMessenger::~XrayFluoMercuryPrimaryGeneratorMessenger()
{
  delete spectrum;
  delete globalCmd;

  G4cout << "XrayFluoMercuryPrimaryGeneratorMessenger deleted" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoMercuryPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 

  if( command == spectrum )
    { XrayFluoAction->SetSpectrum(newValue);} 
  
  if( command == globalCmd )
    { 
      G4bool newGlobalFlag = globalCmd->GetNewBoolValue(newValue);
      XrayFluoAction->SetGlobalFlag(newGlobalFlag);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


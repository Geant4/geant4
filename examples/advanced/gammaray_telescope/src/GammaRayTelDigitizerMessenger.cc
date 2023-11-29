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
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelDigitizerMessenger  ------
//           by F.Longo, G.Santin & R.Giannitrapani (27 nov 2001)
//
// ************************************************************

#include "GammaRayTelDigitizer.hh"
#include "GammaRayTelDigitizerMessenger.hh"

#include "G4SystemOfUnits.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelDigitizerMessenger::GammaRayTelDigitizerMessenger(GammaRayTelDigitizer *GammaRayTelDigit) : GammaRayTelAction(GammaRayTelDigit) {
    constexpr auto ENERGY_DEPOSITION_THRESHOLD_DEFAULT_VALUE{(G4double) 20. * keV};

    thresholdCmd = new G4UIcmdWithADoubleAndUnit("/digitizer/Threshold", this);
    thresholdCmd->SetGuidance("Energy deposition threshold for tracker (TKR) digi generation");
    thresholdCmd->SetParameterName("choice", true);
    thresholdCmd->SetDefaultValue(ENERGY_DEPOSITION_THRESHOLD_DEFAULT_VALUE);
    thresholdCmd->SetRange("Threshold >=0.");
    thresholdCmd->SetUnitCategory("Energy");
    thresholdCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelDigitizerMessenger::~GammaRayTelDigitizerMessenger() {
    delete thresholdCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDigitizerMessenger::SetNewValue(G4UIcommand *command, G4String newValue) {
    if (command == thresholdCmd) {
        GammaRayTelAction->SetThreshold(thresholdCmd->GetNewDoubleValue(newValue));
    }
}

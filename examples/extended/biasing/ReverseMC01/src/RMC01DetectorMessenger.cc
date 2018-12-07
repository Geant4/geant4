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
/// \file biasing/ReverseMC01/src/RMC01DetectorMessenger.cc
/// \brief Implementation of the RMC01DetectorMessenger class
//
//
//////////////////////////////////////////////////////////////
//      Class Name:        RMC01DetectorMessenger
//        Author:               L. Desorgher
//         Organisation:         SpaceIT GmbH
//        Contract:        ESA contract 21435/08/NL/AT
//         Customer:             ESA/ESTEC
//////////////////////////////////////////////////////////////

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RMC01DetectorMessenger.hh"

#include "RMC01DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RMC01DetectorMessenger::RMC01DetectorMessenger(
                                         RMC01DetectorConstruction* GeneralDet)
: G4UImessenger(),
  fTheDetector(GeneralDet),
  fGeneralDir(0),
  fDetDir(0),
  fSetSensitiveVolumeHeightCmd(0),
  fSetSensitiveVolumeRadiusCmd(0),
  fSetShieldingThicknessCmd(0)
{ 
  
  fGeneralDir= new G4UIdirectory("/RMC01/");
  fGeneralDir->SetGuidance(
          "Control of the Geant4 Reverse Monte Carlo example1");
  
  

  fDetDir = new G4UIdirectory("/RMC01/geometry/");
  fDetDir->SetGuidance("Geometry control");
  
  
  fSetSensitiveVolumeHeightCmd = new G4UIcmdWithADoubleAndUnit(
                             "/RMC01/geometry/SetSensitiveVolumeHeight",this);
  fSetSensitiveVolumeHeightCmd->SetGuidance(
                                 "Set the height of the sensitive cylinder");
  fSetSensitiveVolumeHeightCmd->AvailableForStates(G4State_PreInit);
   
  fSetSensitiveVolumeRadiusCmd = new G4UIcmdWithADoubleAndUnit(
                             "/RMC01/geometry/SetSensitiveVolumeRadius",this);
  fSetSensitiveVolumeRadiusCmd->SetGuidance(
                                 "Set the radius of the sensitive cylinder");
  fSetSensitiveVolumeRadiusCmd->AvailableForStates(G4State_PreInit);
  
  fSetShieldingThicknessCmd = new G4UIcmdWithADoubleAndUnit(
                               "/RMC01/geometry/SetShieldingThickness",this);
  fSetShieldingThicknessCmd->SetGuidance(
                        "Set the thickness of the Aluminum Shielding sphere");
  fSetShieldingThicknessCmd->AvailableForStates(G4State_PreInit);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RMC01DetectorMessenger::~RMC01DetectorMessenger()
{
  if (fGeneralDir ) delete fGeneralDir;
  if (fDetDir) delete fDetDir;
  if (fSetSensitiveVolumeHeightCmd) delete  fSetSensitiveVolumeHeightCmd;
  if (fSetSensitiveVolumeRadiusCmd) delete  fSetSensitiveVolumeRadiusCmd;
  if (fSetShieldingThicknessCmd) delete  fSetShieldingThicknessCmd;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RMC01DetectorMessenger::SetNewValue(G4UIcommand* command,
                                                         G4String newValue)
{  if( command == fSetSensitiveVolumeHeightCmd ){ 
          fTheDetector->SetSensitiveVolumeHeight(
                   fSetSensitiveVolumeHeightCmd->GetNewDoubleValue(newValue));
  }
  
  else if( command == fSetSensitiveVolumeRadiusCmd ){ 
          fTheDetector->SetSensitiveVolumeRadius(
                    fSetSensitiveVolumeRadiusCmd->GetNewDoubleValue(newValue));
  }
  
  else if( command == fSetShieldingThicknessCmd ){ 
          fTheDetector->SetShieldingThickness(
                      fSetShieldingThicknessCmd->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

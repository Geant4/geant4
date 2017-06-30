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
// $Id: ExG4DetectorConstruction01Messenger.cc 103117 2017-03-16 14:17:21Z gcosmo $
//
/// \file ExG4DetectorConstruction01Messenger.cc
/// \brief Implementation of the ExG4DetectorConstruction01Messenger class

#include "ExG4DetectorConstruction01Messenger.hh"
#include "ExG4DetectorConstruction01.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADouble.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExG4DetectorConstruction01Messenger::ExG4DetectorConstruction01Messenger(
                           ExG4DetectorConstruction01* detectorConstruction)
 : G4UImessenger(),
   fDetectorConstruction(detectorConstruction),
   fTopDirectory(0),
   fDirectory(0),
   fSetMaterialCmd(0),
   fSetDimensionsCmd(0)
   
{ 
  fTopDirectory = new G4UIdirectory("/ExG4/");
  fTopDirectory->SetGuidance("UI commands of common example classes");
  
  fDirectory = new G4UIdirectory("/ExG4/det/");
  fDirectory->SetGuidance("Detector control");
       
  fSetMaterialCmd 
    = new G4UIcmdWithAString("/ExG4/det/setBoxMaterial",this);
  fSetMaterialCmd->SetGuidance("Select material of the box.");
  fSetMaterialCmd->SetParameterName("BoxMaterial", false);
  fSetMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSetDimensionsCmd 
    = new G4UIcmdWith3VectorAndUnit("/ExG4/det/setBoxDimensions",this);
  fSetDimensionsCmd->SetGuidance("Set box dimensions (in half lentgh).");
  fSetDimensionsCmd->SetParameterName(
                          "BoxDimensionsHx", "BoxDimensionsHy","BoxDimensionsHz", 
                          false);
  //fSetDimensionsCmd->SetUnitCategory("Length");  
  fSetDimensionsCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExG4DetectorConstruction01Messenger::~ExG4DetectorConstruction01Messenger()
{
  delete fTopDirectory;
  delete fDirectory;
  delete fSetMaterialCmd;
  delete fSetDimensionsCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExG4DetectorConstruction01Messenger::SetNewValue(G4UIcommand* command, 
                                                     G4String newValue)
{ 
  if( command == fSetMaterialCmd ) { 
    fDetectorConstruction->SetMaterial(newValue);
  }
  else if( command == fSetDimensionsCmd ) { 
    G4ThreeVector newDimensions 
      = fSetDimensionsCmd->GetNew3VectorValue(newValue);
    fDetectorConstruction
      ->SetDimensions(
          newDimensions.x(), newDimensions.y(), newDimensions.z());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

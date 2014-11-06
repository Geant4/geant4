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
// $Id$
//
/// \file ExG4DetectorConstruction02Messenger.cc
/// \brief Implementation of the ExG4DetectorConstruction02Messenger class

#include "ExG4DetectorConstruction02Messenger.hh"
#include "ExG4DetectorConstruction02.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADouble.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExG4DetectorConstruction02Messenger::ExG4DetectorConstruction02Messenger(
                           ExG4DetectorConstruction02* detectorConstruction)
 : G4UImessenger(),
   fDetectorConstruction(detectorConstruction),
   fTopDirectory(0),
   fDirectory(0),
   fSetBoxMaterialCmd(0),
   fSetWorldMaterialCmd(0),
   fSetBoxDimensionsCmd(0),
   fSetWorldSizeFactorCmd(0)
   
{ 
  fTopDirectory = new G4UIdirectory("/ExG4/");
  fTopDirectory->SetGuidance("UI commands of common example classes");
  
  fDirectory = new G4UIdirectory("/ExG4/det/");
  fDirectory->SetGuidance("Detector control");
       
  fSetBoxMaterialCmd 
    = new G4UIcmdWithAString("/ExG4/det/setBoxMaterial",this);
  fSetBoxMaterialCmd->SetGuidance("Select material of the box.");
  fSetBoxMaterialCmd->SetParameterName("BoxMaterial", false);
  fSetBoxMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
       
  fSetWorldMaterialCmd 
    = new G4UIcmdWithAString("/ExG4/det/setWorldMaterial",this);
  fSetWorldMaterialCmd->SetGuidance("Select material of the world.");
  fSetWorldMaterialCmd->SetParameterName("WorldMaterial", false);
  fSetWorldMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
       
  fSetBoxDimensionsCmd 
    = new G4UIcmdWith3VectorAndUnit("/ExG4/det/setBoxDimensions",this);
  fSetBoxDimensionsCmd->SetGuidance("Set box dimensions (in half lentgh).");
  fSetBoxDimensionsCmd->SetParameterName(
                          "BoxDimensionsHx", "BoxDimensionsHy","BoxDimensionsHz", 
                          false);
  //fSetBoxDimensionsCmd->SetUnitCategory("Length");  
  fSetBoxDimensionsCmd->AvailableForStates(G4State_PreInit);

  fSetWorldSizeFactorCmd 
    = new G4UIcmdWithADouble("/ExG4/det/setWorldSizeFactor",this);
  fSetWorldSizeFactorCmd->SetGuidance(
    "Set the multiplication factor from box dimensions to world dimensions.");
  fSetWorldSizeFactorCmd->SetParameterName("WorldSizeFactor", false);
  fSetWorldSizeFactorCmd->SetRange("WorldSizeFactor >= 1");
  fSetWorldSizeFactorCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExG4DetectorConstruction02Messenger::~ExG4DetectorConstruction02Messenger()
{
  delete fTopDirectory;
  delete fDirectory;
  delete fSetBoxMaterialCmd;
  delete fSetWorldMaterialCmd;
  delete fSetBoxDimensionsCmd;
  delete fSetWorldSizeFactorCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExG4DetectorConstruction02Messenger::SetNewValue(G4UIcommand* command, 
                                                     G4String newValue)
{ 
  if( command == fSetBoxMaterialCmd ) { 
    fDetectorConstruction->SetBoxMaterial(newValue);
  }
  else if( command == fSetWorldMaterialCmd ) { 
    fDetectorConstruction->SetWorldMaterial(newValue);
  }
  else if( command == fSetBoxDimensionsCmd ) { 
    G4ThreeVector newDimensions 
      = fSetBoxDimensionsCmd->GetNew3VectorValue(newValue);
    fDetectorConstruction
      ->SetBoxDimensions(
          newDimensions.x(), newDimensions.y(), newDimensions.z());
  }
  else if ( command == fSetWorldSizeFactorCmd ) { 
    fDetectorConstruction
      ->SetWorldSizeFactor(fSetWorldSizeFactorCmd->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

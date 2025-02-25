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
/// \file DetectorConstructionMessenger.cc
/// \brief Implementation of the DetectorConstructionMessenger class

#include "DetectorConstructionMessenger.hh"

#include "G4UImessenger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstructionMessenger::DetectorConstructionMessenger(DetectorConstruction* det)
: G4UImessenger(), fDetector(det),fTheDetectorDir(nullptr)
{
    fTheDetectorDir = std::make_unique<G4UIdirectory>("/dsbandrepair/det/");
    fTheDetectorDir->SetGuidance("Detector control");

    fWorldDimensionscmd = std::make_unique<G4UIcmdWith3VectorAndUnit>
    ("/dsbandrepair/det/worldBoxSizes",this);
    fWorldDimensionscmd->SetGuidance("Set X, Y, Z sizes for world box");
    fWorldDimensionscmd->SetParameterName("fSizeX","fSizeY","fSizeZ",false);
    fWorldDimensionscmd->AvailableForStates(G4State_PreInit);

    fTheCellDefinitionFilecmd = std::make_unique<G4UIcmdWithAString>
    ("/dsbandrepair/det/celldefinitionfile",this);
    fTheCellDefinitionFilecmd->SetGuidance("Set the path for nucleus definition file");
    fTheCellDefinitionFilecmd->SetParameterName("fTheNucleusDefinitionFile",false);
    fTheCellDefinitionFilecmd->AvailableForStates(G4State_PreInit);

    fTheVoxelDefinitionFilecmd = std::make_unique<G4UIcmdWithAString>
    ("/dsbandrepair/det/voxeldefinitionfile",this);
    fTheVoxelDefinitionFilecmd->SetGuidance("Set the path for voxel type file");
    fTheVoxelDefinitionFilecmd->SetParameterName("fTheVoxelDefinitionFile",false);
    fTheVoxelDefinitionFilecmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstructionMessenger::SetNewValue(G4UIcommand* command, G4String value)
{
    if (command == fTheCellDefinitionFilecmd.get()) {
        fDetector->SetCellDefFilePath(value);
    }
    if (command == fTheVoxelDefinitionFilecmd.get()) {
        fDetector->AddVoxelDefFile(value);
    }
    if (command == fWorldDimensionscmd.get()) {
        auto worldSizesV = fWorldDimensionscmd->GetNew3VectorValue(value);
        fDetector->SetWorldBoxSizes(worldSizesV);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
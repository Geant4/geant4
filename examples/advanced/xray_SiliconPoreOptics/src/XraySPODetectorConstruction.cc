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
/// \file XraySPODetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// Authors: P.Dondero (paolo.dondero@cern.ch), R.Stanzani (ronny.stanzani@cern.ch)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "XraySPODetectorConstruction.hh"
#include "XraySPODetectorMessenger.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"

#include "G4PhysicalVolumeStore.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

XraySPODetectorConstruction::XraySPODetectorConstruction()
: fDetectorMessenger(nullptr)
{
  // materials
  DefineMaterials();
  fReadFile = "";

  // create commands for interactive definition of the calorimeter
  fDetectorMessenger = new XraySPODetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* XraySPODetectorConstruction::Construct()
{
  return ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void XraySPODetectorConstruction::DefineMaterials()
{
  G4NistManager* man = G4NistManager::Instance();

  man->FindOrBuildMaterial("G4_Al");
  man->FindOrBuildMaterial("G4_Si");
  man->FindOrBuildMaterial("G4_Fe");
  man->FindOrBuildMaterial("G4_Cu");
  man->FindOrBuildMaterial("G4_Ge");
  man->FindOrBuildMaterial("G4_Mo");
  man->FindOrBuildMaterial("G4_Ta");
  man->FindOrBuildMaterial("G4_W");
  man->FindOrBuildMaterial("G4_Au");
  man->FindOrBuildMaterial("G4_Pb");
  man->FindOrBuildMaterial("G4_PbWO4");
  man->FindOrBuildMaterial("G4_SODIUM_IODIDE");
  man->FindOrBuildMaterial("G4_AIR");
  man->FindOrBuildMaterial("G4_WATER");

  G4Element* H = man->FindOrBuildElement("H");
  G4Element* O = man->FindOrBuildElement("O");

  G4Material* H2O = new G4Material("Water", 1.000*g/cm3, 2);
  H2O->AddElement(H, 2);
  H2O->AddElement(O, 1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);

  G4double density = universe_mean_density;    //from PhysicalConstants.h
  G4double pressure = 3.e-18*pascal;
  G4double temperature = 2.73*kelvin;
  G4Material* Galactic =  new G4Material("Galactic", 1., 1.008*g/mole, density, kStateGas, temperature, pressure);

  fDefaultMaterial = Galactic;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void XraySPODetectorConstruction::SetReadFile(G4String &input_geometry)
{
  fReadFile = std::move(input_geometry);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* XraySPODetectorConstruction::ConstructDetector()
{
  if (fReadFile == "")
  {
    G4cout << "No geometry selected!" << G4endl;
  }
  else
  {
    auto *InnerRegion = new G4Region("InnerRegion");
    G4cout << "Build geometry from GDML file: " << fReadFile << G4endl;
    G4VPhysicalVolume* fWorldPhysVol;
    fParser.Read(fReadFile);
    fWorldPhysVol = fParser.GetWorldVolume();
    fWorld_log =  fWorldPhysVol->GetLogicalVolume();
    G4int num_volumes = 9;
    auto *inner_volumes = new G4String[num_volumes];
    inner_volumes[0] = "cone_1"; inner_volumes[1] = "cone_2"; inner_volumes[2] = "cone_3"; inner_volumes[3] = "cone_4";
    inner_volumes[4] = "cone_5"; inner_volumes[5] = "cone_6"; inner_volumes[6] = "cone_7"; inner_volumes[7] = "cone_8";

    G4cout << "Adding volumes to the InnerRegion" << G4endl;
    // Add volumes to the InnerRegion
    for (G4int j = 0; j <= num_volumes; j++)
    {
      for (G4int i = 0; i < (G4int)fWorld_log->GetNoDaughters(); i++)
      {
        if (fWorld_log->GetDaughter(i)->GetLogicalVolume()->GetName() == inner_volumes[j])
        {
          InnerRegion->AddRootLogicalVolume(fWorld_log->GetDaughter(i)->GetLogicalVolume());
          G4cout << "Added: " << fWorld_log->GetDaughter(i)->GetLogicalVolume()->GetName() << G4endl;
        }
      }
    }
    fPhysiWorld = fWorldPhysVol;
  }
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


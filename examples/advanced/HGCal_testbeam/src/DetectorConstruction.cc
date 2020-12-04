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
#include "DetectorConstruction.hh"

#include "HGCalTBMaterials.hh"
#include "SiPMSD.hh"
#include "SiliconPixelSD.hh"
#include "DetectorConstruction0.hh"
#include "DetectorConstruction1.hh"
#include "DetectorConstruction2.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4GenericMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4ProductionCuts.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4UserLimits.hh"
#include "G4String.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
    : G4VUserDetectorConstruction(), fConfiguration(-1) {
  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::Construct() {
  // definition of the fMaterials
  fMaterials = new HGCalTBMaterials();
  fMaterials->SetEventDisplayColorScheme();

  /***** Definition of the world = beam line *****/

  // World = Beam line
  G4Box *solidWorld = new G4Box("World", 0.5 * fMaterials->GetBeamLineXY(),
                                0.5 * fMaterials->GetBeamLineXY(),
                                0.5 * fMaterials->GetBeamLineLength());

  G4Material *world_mat = fMaterials->GetAir();
  fLogicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");

  G4VPhysicalVolume *physWorld = new G4PVPlacement(
      0, G4ThreeVector(0., 0., 0.), fLogicWorld, "World", 0, false, 0, true);

  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructHGCal() {

  G4double z0 = -fMaterials->GetBeamLineLength() / 2.;

  std::cout << "Constructing configuration " << fConfiguration << std::endl;

  /*****    START GENERIC PLACEMENT ALGORITHM  FOR THE SETUP  *****/
  for (size_t item_index = 0; item_index < fElementsMap.size(); item_index++) {
    std::string item_type = fElementsMap[item_index].first;
    G4double dz = fElementsMap[item_index].second;
    z0 += dz;

    // places the item at inside the world at z0, z0 is incremented by the
    // item's thickness
    fMaterials->PlaceItemInLogicalVolume(item_type, z0, fLogicWorld);
  }

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4UImanager *UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand("/vis/drawVolume");
  UImanager->ApplyCommand("/vis/viewer/set/targetPoint 0 0 " +
                          std::to_string(fVisViewpoint / CLHEP::m) + " m");
  UImanager->ApplyCommand("/vis/scene/add/trajectories smooth");
  UImanager->ApplyCommand("/vis/scene/add/hits");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField() {
  G4SDManager *sdman = G4SDManager::GetSDMpointer();

  SiliconPixelSD *sensitiveSilicon = new SiliconPixelSD(
      (fMaterials->GetSiPixelLogical()->GetName() + "_sensitive").c_str());
  sdman->AddNewDetector(sensitiveSilicon);
  fMaterials->GetSiPixelLogical()->SetSensitiveDetector(sensitiveSilicon);

  SiPMSD *sensitiveSiPM = new SiPMSD(
      (fMaterials->GetAHCALSiPMlogical()->GetName() + "_sensitive").c_str());
  sdman->AddNewDetector(sensitiveSiPM);
  fMaterials->GetAHCALSiPMlogical()->SetSensitiveDetector(sensitiveSiPM);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SelectConfiguration(G4int val) {

  if (fConfiguration != -1) {
    G4ExceptionDescription msg;
    msg << "Configuration " << fConfiguration << " is already placed.\n"
        << "Configuration can be set only once. Please restart (and\n"
        << "edit your macro if necessary).\n";
    G4Exception("DetectorConstruction::SelectConfiguration()", "MultipleConfig",
                JustWarning, msg);
    return;
  }

  fVisViewpoint = 0;
  if (val == 0)
    DetectorConstruction0(fElementsMap, fVisViewpoint);
  else if (val == 1)
    DetectorConstruction1(fElementsMap, fVisViewpoint);
  else if (val == 2)
    DetectorConstruction2(fElementsMap, fVisViewpoint);
  else {
    G4ExceptionDescription msg;
    msg << "Configuration " << val << " is not implemented.\n"
        << "Choose between configuration 0, 1, and 2.\n";
    G4Exception("DetectorConstruction::SelectConfiguration()", "WrongConfig",
                JustWarning, msg);
    return;
  }
  fConfiguration = val;

  ConstructHGCal();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetStepSizeSilicon(G4double val) {
  // setting the step size in silicon:
  G4double maxTrackLength = val * 0.001 * CLHEP::mm;
  fMaterials->GetSiPixelLogical()->SetUserLimits(
      new G4UserLimits(0, maxTrackLength));

  G4Region *reg = fMaterials->GetSiPixelLogical()->GetRegion();
  G4ProductionCuts *cuts = new G4ProductionCuts;
  cuts->SetProductionCut(maxTrackLength);
  reg->SetProductionCuts(cuts);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineCommands() {
  // define command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this, "/HGCalTestbeam/setup/",
                                      "Configuration specifications");

  // configuration command
  auto &configCmd = fMessenger->DeclareMethod(
      "configuration", &DetectorConstruction::SelectConfiguration,
      "Select the configuration (0 for HGCal test beam, 1 for same HGCal"
      " with beamline (upstream material), or 2 for simple test setup)");
  configCmd.SetParameterName("index", true);
  configCmd.SetDefaultValue("0");

  auto &SiStepSizeCmd = fMessenger->DeclareMethod(
      "stepSilicon", &DetectorConstruction::SetStepSizeSilicon,
      "Maximum step size in silicon pixels, unit: microns");
  SiStepSizeCmd.SetParameterName("size", true);
  SiStepSizeCmd.SetDefaultValue("30.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

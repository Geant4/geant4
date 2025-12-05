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
/// \file GB03DetectorConstruction.cc
/// \brief Implementation of the GB03DetectorConstruction class

#include "GB03DetectorConstruction.hh"

#include "GB03BOptrGeometryBasedBiasing.hh"
#include "GB03HodoSD.hh"

#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4GenericMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

G4int GB03DetectorConstruction::fNumberOfLayers = 40;
G4ThreadLocal G4bool GB03DetectorConstruction::fConstructedSDandField = false;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB03DetectorConstruction::GB03DetectorConstruction(G4bool bf)
  : G4VUserDetectorConstruction(),
    fTotalThickness(2.0 * m),
    fLayerThickness(0.),
    fConstructed(false),
    fWorldMaterial(nullptr),
    fAbsorberMaterial(nullptr),
    fGapMaterial(nullptr),
    fLayerSolid(nullptr),
    fGapSolid(nullptr),
    fCalorLogical(nullptr),
    fLayerLogical(nullptr),
    fGapLogical(nullptr),
    fWorldPhysical(nullptr),
    fLayerPhysical(nullptr),
    fGapPhysical(nullptr),
    fDetectorMessenger(nullptr),
    fVerboseLevel(1),
    fBiasingFlag(bf)
{
  fLayerThickness = fTotalThickness / fNumberOfLayers;
  fCalName = "Calor";

  DefineMaterials();

  // Material List
  G4String matList;
  const G4MaterialTable* matTbl = G4Material::GetMaterialTable();
  for (size_t i = 0; i < G4Material::GetNumberOfMaterials(); i++) {
    matList += (*matTbl)[i]->GetName();
    matList += " ";
  }

  // Detector Messenger
  fDetectorMessenger = new G4GenericMessenger(this, "/GB03/", "UI commands of this example");
  fDetectorMessenger
    ->DeclareMethod("setAbsMat", &GB03DetectorConstruction::SetAbsorberMaterial,
                    "Select Material of the Absorber.")
    .SetParameterName("material", false)
    .SetStates(G4State_Idle)
    .SetCandidates(matList);

  fDetectorMessenger
    ->DeclareMethod("setGapMat", &GB03DetectorConstruction::SetGapMaterial,
                    "Select Material of the Gap.")
    .SetParameterName("choice", false)
    .SetStates(G4State_Idle)
    .SetCandidates(matList);

  fDetectorMessenger
    ->DeclareMethod("numberOfLayers", &GB03DetectorConstruction::SetNumberOfLayers,
                    "Set number of layers.")
    .SetParameterName("nl", false)
    .SetStates(G4State_Idle)
    .SetRange("nl>0");

  fDetectorMessenger
    ->DeclareMethod("verbose", &GB03DetectorConstruction::SetVerboseLevel, "Set verbosity level.")
    .SetParameterName("verbose", false)
    .SetStates(G4State_Idle)
    .SetRange("verbose>=0");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB03DetectorConstruction::~GB03DetectorConstruction()
{
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* GB03DetectorConstruction::Construct()
{
  if (!fConstructed) {
    fConstructed = true;
    SetupGeometry();
  }
  if (GetVerboseLevel() > 0) {
    PrintCalorParameters();
  }

  return fWorldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB03DetectorConstruction::ConstructSDandField()
{
  if (!fConstructedSDandField) {
    fConstructedSDandField = true;
    SetupDetectors();
    if (fBiasingFlag) {
      SetupBiasing();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB03DetectorConstruction::DefineMaterials()
{
  G4String name, symbol;  // a=mass of a mole;
  G4double a, z, density;  // z=mean number of protons;
  G4int iz;  // iz=number of protons  in an isotope;
  G4int n;  // n=number of nucleons in an isotope;

  G4int ncomponents, natoms;
  G4double abundance, fractionmass;
  G4double temperature, pressure;

  //
  // define Elements
  //

  a = 1.01 * g / mole;
  auto H = new G4Element(name = "Hydrogen", symbol = "H", z = 1., a);

  a = 12.01 * g / mole;
  auto C = new G4Element(name = "Carbon", symbol = "C", z = 6., a);

  a = 14.01 * g / mole;
  auto N = new G4Element(name = "Nitrogen", symbol = "N", z = 7., a);

  a = 16.00 * g / mole;
  auto O = new G4Element(name = "Oxygen", symbol = "O", z = 8., a);

  //
  // define an Element from isotopes, by relative abundance
  //

  auto U5 = new G4Isotope(name = "U235", iz = 92, n = 235, a = 235.01 * g / mole);
  auto U8 = new G4Isotope(name = "U238", iz = 92, n = 238, a = 238.03 * g / mole);

  auto U = new G4Element(name = "enriched Uranium", symbol = "U", ncomponents = 2);
  U->AddIsotope(U5, abundance = 90. * perCent);
  U->AddIsotope(U8, abundance = 10. * perCent);

  //
  // define simple materials
  //

  new G4Material(name = "Aluminium", z = 13., a = 26.98 * g / mole, density = 2.700 * g / cm3);
  new G4Material(name = "Silicon", z = 14., a = 28.09 * g / mole, density = 2.33 * g / cm3);
  new G4Material(name = "Iron", z = 26., a = 55.85 * g / mole, density = 7.87 * g / cm3);
  new G4Material(name = "ArgonGas", z = 18., a = 39.95 * g / mole, density = 1.782 * mg / cm3);
  new G4Material(name = "He", z = 2., a = 4.0 * g / mole, density = 0.1786e-03 * g / cm3);

  density = 1.390 * g / cm3;
  a = 39.95 * g / mole;
  new G4Material(name = "liquidArgon", z = 18., a, density);

  density = 11.35 * g / cm3;
  a = 207.19 * g / mole;
  auto Pb = new G4Material(name = "Lead", z = 82., a, density);

  //
  // define a material from elements.   case 1: chemical molecule
  //

  density = 1.000 * g / cm3;
  auto H2O = new G4Material(name = "Water", density, ncomponents = 2);
  H2O->AddElement(H, natoms = 2);
  H2O->AddElement(O, natoms = 1);

  density = 1.032 * g / cm3;
  auto Sci = new G4Material(name = "Scintillator", density, ncomponents = 2);
  Sci->AddElement(C, natoms = 9);
  Sci->AddElement(H, natoms = 10);

  //
  // define a material from elements.   case 2: mixture by fractional mass
  //

  density = 1.290 * mg / cm3;
  auto Air = new G4Material(name = "Air", density, ncomponents = 2);
  Air->AddElement(N, fractionmass = 0.7);
  Air->AddElement(O, fractionmass = 0.3);

  //
  // examples of vacuum
  //

  density = universe_mean_density;
  pressure = 3.e-18 * pascal;
  temperature = 2.73 * kelvin;
  auto Vacuum = new G4Material(name = "Galactic", z = 1., a = 1.01 * g / mole, density, kStateGas,
                               temperature, pressure);

  if (GetVerboseLevel() > 1) {
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  }

  // default materials of the calorimeter
  fWorldMaterial = Vacuum;
  fAbsorberMaterial = Pb;
  fGapMaterial = Sci;
  fHodoMaterial = Sci;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB03DetectorConstruction::SetupGeometry()
{
  //
  // World
  //
  G4VSolid* worldSolid = new G4Box("World", 2. * m, 2. * m, fTotalThickness * 2.);
  auto* WorldLogical = new G4LogicalVolume(worldSolid, fWorldMaterial, "World");
  fWorldPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(), WorldLogical, "World", nullptr, false, 0);

  //
  // Calorimeter
  //
  G4VSolid* calorSolid = new G4Box("Calor", 0.5 * m, 0.5 * m, fTotalThickness / 2.);
  fCalorLogical = new G4LogicalVolume(calorSolid, fAbsorberMaterial, fCalName);
  new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), fCalorLogical, fCalName, WorldLogical,
                    false, 0);

  //
  // Layers --- as absorbers
  //
  fLayerSolid = new G4Box("Layer", 0.5 * m, 0.5 * m, fLayerThickness / 2.);
  fLayerLogical = new G4LogicalVolume(fLayerSolid, fAbsorberMaterial, fCalName + "_LayerLog");
  fLayerPhysical = new G4PVReplica(fCalName + "_Layer", fLayerLogical, fCalorLogical, kZAxis,
                                   fNumberOfLayers, fLayerThickness);

  //
  // Gap
  //
  fGapSolid = new G4Box("Gap", 0.5 * m, 0.5 * m, fLayerThickness / 4.);
  fGapLogical = new G4LogicalVolume(fGapSolid, fGapMaterial, fCalName + "_Gap");
  fGapPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0., 0., fLayerThickness / 4.),
                                   fGapLogical, fCalName + "_gap", fLayerLogical, false, 0);

  //
  // Hodoscope plane
  //
  G4VSolid* measSolid = new G4Box("Hodo", 0.5 * m, 0.5 * m, 5 * cm);
  fHodoLogical = new G4LogicalVolume(measSolid, fHodoMaterial, "Hodo_l");
  new G4PVPlacement(nullptr, G4ThreeVector(0., 0., fTotalThickness / 2. + 5 * cm), fHodoLogical,
                    "Hodo_p", WorldLogical, false, 0);
  //
  // Visualization attributes
  //
  WorldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto simpleBoxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
  auto GapBoxVisAtt = new G4VisAttributes(G4Colour::Blue());
  // auto LayerBoxVisAtt = new G4VisAttributes(G4Colour::Green());
  auto HodoBoxVisAtt = new G4VisAttributes(G4Colour::Red());
  simpleBoxVisAtt->SetVisibility(true);
  fCalorLogical->SetVisAttributes(simpleBoxVisAtt);
  fLayerLogical->SetVisAttributes(simpleBoxVisAtt);
  fGapLogical->SetVisAttributes(GapBoxVisAtt);
  fHodoLogical->SetVisAttributes(HodoBoxVisAtt);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB03DetectorConstruction::SetupDetectors()
{
  // ---------------------------------------------------------------------------------
  // -- Attach sensitive detector to print information on particle exiting the shield:
  // ---------------------------------------------------------------------------------
  // -- create and register sensitive detector code module:
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4VSensitiveDetector* sd = new GB03HodoSD("HodoSD");
  SDman->AddNewDetector(sd);
  // -- Fetch volume for sensitivity and attach sensitive module to it:
  fHodoLogical->SetSensitiveDetector(sd);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB03DetectorConstruction::SetupBiasing()
{
  auto biasingOperator = new GB03BOptrGeometryBasedBiasing();
  biasingOperator->AttachTo(fLayerLogical);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB03DetectorConstruction::PrintCalorParameters() const
{
  G4cout << "--------------------------------------------------------" << G4endl;
  G4cout << " Absorber is made of " << fAbsorberMaterial->GetName() << G4endl << " Gap is made of "
         << fGapMaterial->GetName() << G4endl
         << "--------------------------------------------------------" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB03DetectorConstruction::SetAbsorberMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) {
    fAbsorberMaterial = pttoMaterial;
    if (fConstructed) {
      fCalorLogical->SetMaterial(fAbsorberMaterial);
      fLayerLogical->SetMaterial(fAbsorberMaterial);
    }
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    if (GetVerboseLevel() > 1) {
      PrintCalorParameters();
    }
  }
  else {
    G4cerr << materialChoice << " is not defined. - Command is ignored." << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String GB03DetectorConstruction::GetAbsorberMaterial() const
{
  return fAbsorberMaterial->GetName();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB03DetectorConstruction::SetGapMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) {
    fGapMaterial = pttoMaterial;
    if (fConstructed) {
      fGapLogical->SetMaterial(fGapMaterial);
    }
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    if (GetVerboseLevel() > 1) {
      PrintCalorParameters();
    }
  }
  else {
    G4cerr << materialChoice << " is not defined. - Command is ignored." << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String GB03DetectorConstruction::GetGapMaterial() const
{
  return fGapMaterial->GetName();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB03DetectorConstruction::SetNumberOfLayers(G4int nl)
{
  fNumberOfLayers = nl;
  fLayerThickness = fTotalThickness / fNumberOfLayers;
  if (!fConstructed) return;

  fLayerSolid->SetZHalfLength(fLayerThickness / 2.);
  fGapSolid->SetZHalfLength(fLayerThickness / 4.);

  fCalorLogical->RemoveDaughter(fLayerPhysical);
  delete fLayerPhysical;
  fLayerPhysical = new G4PVReplica(fCalName + "_Layer", fLayerLogical, fCalorLogical, kZAxis,
                                   fNumberOfLayers, fLayerThickness);
  fGapPhysical->SetTranslation(G4ThreeVector(0., 0., fLayerThickness / 4.));

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

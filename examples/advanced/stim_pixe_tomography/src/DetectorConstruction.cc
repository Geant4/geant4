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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"

#include "G4AutoDelete.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Ellipsoid.hh"
#include "G4GeometryManager.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4UniformMagField.hh"
#include "G4UnitsTable.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"

#include "DetectorMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(),
    fAbsorberMaterial(nullptr),
    fWorldMaterial(nullptr),
    fSolidWorld(nullptr),
    fLogicWorld(nullptr),
    fPhysiWorld(nullptr),
    fSolidAbsorber(nullptr),
    fLogicAbsorber(nullptr),
    fPhysiAbsorber(nullptr),
    material1(nullptr),
    material2(nullptr),
    material3(nullptr),
    material4(nullptr),
    material5(nullptr),
    material6(nullptr),
    material_GDP(nullptr),
    ellipse1(nullptr),
    logicEllipse1(nullptr),
    physiEllipse1(nullptr),
    ellipse2(nullptr),
    logicEllipse2(nullptr),
    physiEllipse2(nullptr),
    ellipse3(nullptr),
    logicEllipse3(nullptr),
    physiEllipse3(nullptr),
    ellipse4(nullptr),
    logicEllipse4(nullptr),
    physiEllipse4(nullptr),
    ellipse5(nullptr),
    logicEllipse5(nullptr),
    physiEllipse5(nullptr),
    ellipse6(nullptr),
    logicEllipse6(nullptr),
    physiEllipse6(nullptr),
    solid_GDP(nullptr),
    logic_GDP(nullptr),
    physi_GDP(nullptr),
    fDetectorMessenger(nullptr)
{
  phantom_type = 1;

  DefineMaterials();

  if (phantom_type == 1) {
    // default parameter values of the box
    fAbsorberThickness = 40. * um;
    fAbsorberSizeYZ = 40. * um;

    //    SetAbsorberMaterial("G4_Cu");
    SetAbsorberMaterial("Body_10times");
    //    SetAbsorberMaterial("Body_real");
  }
  fXposAbs = 0. * cm;
  ComputeGeomParameters();

  // materials
  SetWorldMaterial("G4_Galactic");

  // create commands for interactive definition of the box
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  //
  // define Elements
  //
  auto H = new G4Element("Hydrogen", "H", 1, 1.01 * g / mole);
  auto C = new G4Element("Carbon", "C", 6, 12.01 * g / mole);
  auto N = new G4Element("Nitrogen", "N", 7, 14.01 * g / mole);
  auto O = new G4Element("Oxygen", "O", 8, 16.00 * g / mole);
  auto P = new G4Element("Phosphorus", "P", 15, 30.97 * g / mole);
  auto S = new G4Element("Sulphur", "S", 16, 32.07 * g / mole);
  auto Cl = new G4Element("Chlorine", "Cl", 17, 35.45 * g / mole);
  auto K = new G4Element("Potassium", "K", 19, 39.10 * g / mole);
  auto Ca = new G4Element("Calcium", "Ca", 20, 40.08 * g / mole);
  auto Ti = new G4Element("Titanium", "Ti", 22, 47.87 * g / mole);
  auto Ge = new G4Element("Germanium", "Ge", 32., 72.59 * g / mole);

  //
  // define biological materials
  //
  auto BioMaterial = new G4Material("C10H17O3N2", 1.218 * g / cm3, 4);
  BioMaterial->AddElement(C, 10);
  BioMaterial->AddElement(H, 17);
  BioMaterial->AddElement(O, 3);
  BioMaterial->AddElement(N, 2);

  // FIRST ELLIPSOID	: external surface of the worm ("skin")
  // Original (real) composition in dry mass (the worm is dehydrated by
  // lyophilization) Each element content is defined by its mass fraction (must
  // be normalised to 1)
  auto Skin_real = new G4Material("Skin_real", 0.49729 * g / cm3, 5);
  Skin_real->AddElement(P, 0.47 * perCent);
  Skin_real->AddElement(S, 0.21 * perCent);
  Skin_real->AddElement(Cl, 0.05 * perCent);
  Skin_real->AddElement(K, 0.48 * perCent);
  Skin_real->AddMaterial(BioMaterial, 98.79 * perCent);

  // Skin1: mass fractions of mineral elements 1O times more than original
  // (real) mass fractions
  auto Skin_10times = new G4Material("Skin_10times", 0.49729 * g / cm3, 5);
  Skin_10times->AddElement(P, 4.71 * perCent);
  Skin_10times->AddElement(S, 2.10 * perCent);
  Skin_10times->AddElement(Cl, 0.46 * perCent);
  Skin_10times->AddElement(K, 4.77 * perCent);
  Skin_10times->AddMaterial(BioMaterial, 87.96 * perCent);

  // SECOND ELLIPSOID	: "body" of the worm, limited by the internal surface of
  // the "skin" Original (real) composition in dry mass (the worm is dehydrated
  // by lyophilization) Each element content is defined by its mass fraction
  // (must be normalised to 1)
  auto Body_real = new G4Material("Body_real", 0.40853 * g / cm3, 2);
  Body_real->AddElement(P, 0.72 * perCent);
  Body_real->AddMaterial(BioMaterial, 99.28 * perCent);

  // Body_10times: mass fractions of mineral elements 1O times more than
  // original (real) mass fractions
  auto Body_10times = new G4Material("Body_10times", 0.40853 * g / cm3, 2);
  Body_10times->AddElement(P, 7.23 * perCent);
  Body_10times->AddMaterial(BioMaterial, 92.77 * perCent);

  // THIRD ELLIPSOID : Nucleus of cell 1 (contained inside ellipsoid 2)
  // Original (real) composition in dry mass (the worm is dehydrated by
  // lyophilization) Each element content is defined by its mass fraction (must
  // be normalised to 1)
  auto Cell1_real = new G4Material("Cell1_real", 0.66262 * g / cm3, 3);
  Cell1_real->AddElement(P, 0.57 * perCent);
  Cell1_real->AddElement(Ca, 1.55 * perCent);
  Cell1_real->AddMaterial(BioMaterial, 97.88 * perCent);

  // Cell1_10times: mass fractions of mineral elements 1O times more than
  // original (real) mass fractions
  auto Cell1_10times = new G4Material("Cell1_10times", 0.66262 * g / cm3, 3);
  Cell1_10times->AddElement(P, 5.66 * perCent);
  Cell1_10times->AddElement(Ca, 15.49 * perCent);
  Cell1_10times->AddMaterial(BioMaterial, 78.85 * perCent);

  // FOURTH ELLIPSOID : Nucleus of cell 2 (contained inside ellipsoid 2)
  // Original (real) composition in dry mass (the worm is dehydrated by
  // lyophilization) Each element content is defined by its mass fraction (must
  // be normalised to 1)
  auto Cell2_real = new G4Material("Cell2_real", 0.60227 * g / cm3, 3);
  Cell2_real->AddElement(P, 0.61 * perCent);
  Cell2_real->AddElement(Ca, 1.44 * perCent);
  Cell2_real->AddMaterial(BioMaterial, 97.95 * perCent);

  // Cell2_10times: mass fractions of mineral elements 1O times more than
  // original (real) mass fractions
  auto Cell2_10times = new G4Material("Cell2_10times", 0.60227 * g / cm3, 3);
  Cell2_10times->AddElement(P, 6.10 * perCent);
  Cell2_10times->AddElement(Ca, 14.45 * perCent);
  Cell2_10times->AddMaterial(BioMaterial, 79.45 * perCent);

  // FIFTH ELLIPSOID: intestine of the worm
  // Original (real) composition in dry mass (the worm is dehydrated by
  // lyophilization) Each element content is defined by its mass fraction (must
  // be normalised to 1)
  auto Intestine_real = new G4Material("Intestine_real", 0.54058 * g / cm3, 4);
  Intestine_real->AddElement(S, 0.12 * perCent);
  Intestine_real->AddElement(Cl, 0.02 * perCent);
  Intestine_real->AddElement(K, 0.20 * perCent);
  Intestine_real->AddMaterial(BioMaterial, 99.66 * perCent);

  // Intestine_10times: mass fractions of mineral elements 1O times more than
  // original (real) mass fractions
  auto Intestine_10times = new G4Material("Intestine_10times", 0.54058 * g / cm3, 4);
  Intestine_10times->AddElement(S, 1.17 * perCent);
  Intestine_10times->AddElement(Cl, 0.23 * perCent);
  Intestine_10times->AddElement(K, 1.99 * perCent);
  Intestine_10times->AddMaterial(BioMaterial, 96.61 * perCent);

  // SIXTH ELLIPSOID: Nucleus of cell Ti (contained inside ellipsoid 5)
  // Original (real) composition in dry mass (the worm is dehydrated by
  // lyophilization) Each element content is defined by its mass fraction (must
  // be normalised to 1)
  auto CellTi_real = new G4Material("CellTi_real", 0.75138 * g / cm3, 5);
  CellTi_real->AddElement(S, 0.11 * perCent);
  CellTi_real->AddElement(Cl, 0.02 * perCent);
  CellTi_real->AddElement(K, 0.18 * perCent);
  CellTi_real->AddElement(Ti, 0.05 * perCent);
  CellTi_real->AddMaterial(BioMaterial, 99.64 * perCent);

  // CellTi_200times: mass fractions of mineral elements 1O times more than
  // original (real) mass fractions except for Ti which is multiplied by 200
  // times
  auto CellTi_200times = new G4Material("CellTi_200times", 0.75138 * g / cm3, 5);
  CellTi_200times->AddElement(S, 1.05 * perCent);
  CellTi_200times->AddElement(Cl, 0.22 * perCent);
  CellTi_200times->AddElement(K, 1.79 * perCent);
  CellTi_200times->AddElement(Ti, 10.89 * perCent);  // Ti density 200 times
  CellTi_200times->AddMaterial(BioMaterial, 86.05 * perCent);

  // CellTi_1000times: mass fractions of mineral elements 10 times more than
  // original (real) mass fractions except for Ti which is multiplied by 1000
  // times
  auto CellTi_1000times = new G4Material("CellTi_1000times", 0.75138 * g / cm3, 5);
  CellTi_1000times->AddElement(S, 1.05 * perCent);
  CellTi_1000times->AddElement(Cl, 0.22 * perCent);
  CellTi_1000times->AddElement(K, 1.79 * perCent);
  CellTi_1000times->AddElement(Ti, 54.47 * perCent);  // Ti density 1000 times
  CellTi_1000times->AddMaterial(BioMaterial, 42.47 * perCent);

  //
  // define simple materials  Cl-doped polystyrene capsules (PS)
  //
  // Polystyrene
  auto Cl_Polystyrene = new G4Material("Cl-doped Polystyrene", 1.18425 * g / cm3, 3);
  Cl_Polystyrene->AddElement(C, 97);
  Cl_Polystyrene->AddElement(H, 97);
  Cl_Polystyrene->AddElement(Cl, 6);

  // Ge-doped glow discharge polymer (GDP)
  auto GDP = new G4Material("GDP", 1.08 * g / cm3, 3);
  GDP->AddElement(C, 76430);
  GDP->AddElement(H, 107002);
  GDP->AddElement(Ge, 744);

  //
  // define a material from elements.
  //
  auto H2O = new G4Material("Water", 1.000 * g / cm3, 2);
  H2O->AddElement(H, 2);
  H2O->AddElement(O, 1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78 * eV);

  auto Air = new G4Material("Air", 1.290 * mg / cm3, 2);
  Air->AddElement(N, 0.7);
  Air->AddElement(O, 0.3);

  //
  // example of vacuum
  //
  G4double density = universe_mean_density;  // from PhysicalConstants.h
  G4double pressure = 3.e-18 * pascal;
  G4double temperature = 2.73 * kelvin;
  new G4Material("Galactic", 1, 1.01 * g / mole, density, kStateGas, temperature, pressure);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ComputeGeomParameters()
{
  if (phantom_type == 1) {
    // Compute derived parameters of the box
    fXstartAbs = fXposAbs - 0.5 * fAbsorberThickness;
    fXendAbs = fXposAbs + 0.5 * fAbsorberThickness;

    G4double xmax = std::max(std::abs(fXstartAbs), std::abs(fXendAbs));
    fWorldSizeX = 4 * xmax;
    fWorldSizeYZ = 2 * fAbsorberSizeYZ;

    if (nullptr != fPhysiWorld) {
      ChangeGeometry();
    }
  }
  else if (phantom_type == 2) {
    fWorldSizeX = 1 * mm;
    fWorldSizeYZ = 1 * mm;
  }
  else if (phantom_type == 3) {
    fWorldSizeX = 1.5 * mm;
    fWorldSizeYZ = 1.5 * mm;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  if (nullptr != fPhysiWorld) {
    return fPhysiWorld;
  }
  // World
  //
  fSolidWorld = new G4Box("World",  // its name
                          fWorldSizeX / 2, fWorldSizeYZ / 2, fWorldSizeYZ / 2);  // its size

  fLogicWorld = new G4LogicalVolume(fSolidWorld,  // its solid
                                    fWorldMaterial,  // its material
                                    "World");  // its name

  fPhysiWorld = new G4PVPlacement(nullptr,  // no rotation
                                  G4ThreeVector(0., 0., 0.),  // at (0,0,0)
                                  fLogicWorld,  // its logical volume
                                  "World",  // its name
                                  nullptr,  // its mother  volume
                                  false,  // no boolean operation
                                  0);  // copy number
  if (phantom_type == 1) {
    Construct_Phantom1();
  }
  else if (phantom_type == 2) {
    Construct_Phantom2();
  }
  else if (phantom_type == 3) {
    Construct_Phantom3();
  }

  //    fLogicWorld->SetVisAttributes (G4VisAttributes::GetInvisible());

  // always return the physical World
  //
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintGeomParameters()
{
  if (phantom_type == 1) {
    G4cout << "\n" << fWorldMaterial << G4endl;
    G4cout << "\n" << fAbsorberMaterial << G4endl;

    G4cout << "\n The  WORLD   is made of " << G4BestUnit(fWorldSizeX, "Length") << " of "
           << fWorldMaterial->GetName();
    G4cout << ". The transverse size (YZ) of the world is " << G4BestUnit(fWorldSizeYZ, "Length")
           << G4endl;
    G4cout << " The ABSORBER is made of " << G4BestUnit(fAbsorberThickness, "Length") << " of "
           << fAbsorberMaterial->GetName();
    G4cout << ". The transverse size (YZ) is " << G4BestUnit(fAbsorberSizeYZ, "Length") << G4endl;
    G4cout << " X position of the middle of the absorber " << G4BestUnit(fXposAbs, "Length");
    G4cout << G4endl;
  }
  else if (phantom_type == 2) {
    G4cout << "The upper part of C.elegans is made of 6 ellipsoids" << G4endl;
  }
  else if (phantom_type == 3) {
    G4cout << "A microsphere inertial confinement fusion target is contructed, made of: " << G4endl;
    G4cout << "\n" << material_GDP << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberMaterial(const G4String& materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fAbsorberMaterial != pttoMaterial) {
    fAbsorberMaterial = pttoMaterial;
    if (fLogicAbsorber) {
      fLogicAbsorber->SetMaterial(fAbsorberMaterial);
    }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(const G4String& materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fWorldMaterial != pttoMaterial) {
    fWorldMaterial = pttoMaterial;
    if (fLogicWorld) {
      fLogicWorld->SetMaterial(fWorldMaterial);
    }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberThickness(G4double val)
{
  fAbsorberThickness = val;
  ComputeGeomParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberSizeYZ(G4double val)
{
  fAbsorberSizeYZ = val;
  ComputeGeomParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldSizeX(G4double val)
{
  fWorldSizeX = val;
  ComputeGeomParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldSizeYZ(G4double val)
{
  fWorldSizeYZ = val;
  ComputeGeomParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberXpos(G4double val)
{
  fXposAbs = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::ConstructSDandField()
{
  if (fFieldMessenger.Get() == nullptr) {
    // Create global magnetic field messenger.
    // Uniform magnetic field is then created automatically if
    // the field value is not zero.
    G4ThreeVector fieldValue = G4ThreeVector();
    auto msg = new G4GlobalMagFieldMessenger(fieldValue);
    // msg->SetVerboseLevel(1);
    G4AutoDelete::Register(msg);
    fFieldMessenger.Put(msg);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ChangeGeometry()
{
  if (phantom_type == 1) {
    fSolidWorld->SetXHalfLength(fWorldSizeX * 0.5);
    fSolidWorld->SetYHalfLength(fWorldSizeYZ * 0.5);
    fSolidWorld->SetZHalfLength(fWorldSizeYZ * 0.5);

    fSolidAbsorber->SetXHalfLength(fAbsorberThickness * 0.5);
    fSolidAbsorber->SetYHalfLength(fAbsorberSizeYZ * 0.5);
    fSolidAbsorber->SetZHalfLength(fAbsorberSizeYZ * 0.5);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::SetPhantomType(G4int value)
{
  phantom_type = value;
  ComputeGeomParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::Construct_Phantom1()
{
  fSolidAbsorber =
    new G4Box("Absorber", fAbsorberThickness / 2, fAbsorberSizeYZ / 2, fAbsorberSizeYZ / 2);

  fLogicAbsorber = new G4LogicalVolume(fSolidAbsorber,  // its solid
                                       fAbsorberMaterial,  // its material
                                       "Absorber");  // its name

  G4RotationMatrix rm;
  G4double angle = 0 * deg;
  rm.rotateZ(angle);

  fPhysiAbsorber = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(fXposAbs * um, 0., 0.)),
                                     fLogicAbsorber,  // its logical volume
                                     "Absorber",  // its name
                                     fLogicWorld,  // its mother
                                     false,  // no boolean operation
                                     0);  // copy number

  auto logVisAtt = new G4VisAttributes(G4Colour(0.670588, 0.333333, 0.862745));
  logVisAtt->SetForceSolid(true);
  fLogicAbsorber->SetVisAttributes(logVisAtt);

  PrintGeomParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::Construct_Phantom2()
{
  G4NistManager* nist = G4NistManager::Instance();
  // PHANTOM 3: C. ELEGANS WORM **********************************
  // made of 6 ellipsoids

  // FIRST ELLIPSOID	: external surface of the worm ("skin")

  material1 = nist->FindOrBuildMaterial("Skin_10times");
  G4ThreeVector center1 =
    G4ThreeVector(fXposAbs, 0,
                  0);  // position of the center of the first ellipsoid is (fXposAbs,0,0)
  // if we change fXposAbs, all the phantom will be translated along X axis
  // as the position of the center of all other ellipsoids is defined according
  // to fXposAbs
  G4double halfx1 = 20.61 * um;  // Semiaxis in X
  G4double halfy1 = 21.42 * um;  // Semiaxis in Y
  G4double halfz1 = 187.82 * um;  // Semiaxis in Z
  G4double pzTopCut1 = 187.82 * um;  // here the top of the ellipsoid is not cut away
  G4double pzBottomCut1 = 0 * um;  // here the bottom of the ellipsoid is cut
                                   // away below z = 0 (original configuration)
  // G4double pzBottomCut1 = 18.07*um;  // test with position of the slice of
  // interest for vis.mac visualisation of the slice z = 18.07 µm

  ellipse1 = new G4Ellipsoid("Ellipse1", halfx1, halfy1, halfz1, pzBottomCut1,
                             pzTopCut1);  // cut ellipsoid (original configuration)

  logicEllipse1 = new G4LogicalVolume(ellipse1,  // its solid
                                      material1,  // its material
                                      "Ellipse1");  // its name

  physiEllipse1 = new G4PVPlacement(nullptr,  // no rotation
                                    center1,  // its position
                                    logicEllipse1,  // its logical volume
                                    "Ellipse1",  // its name
                                    fLogicWorld,  // its mother
                                    false,  // no boolean operation
                                    0);  // copy number

  // SECOND ELLIPSOID	: "body" of the worm, limited by the internal surface of
  // the "skin"
  material2 = nist->FindOrBuildMaterial("Body_10times");
  G4ThreeVector center2 = G4ThreeVector(-0.39 * um, 0, 0);  // position according to ellipsoid 1
  G4double halfx2 = 18.61 * um;  // Semiaxis in X
  G4double halfy2 = 19.01 * um;  //  Semiaxis in Y
  G4double halfz2 = 186.64 * um;  // Semiaxis in Z

  G4double pzTopCut2 = 186.64 * um;  // here the top of the ellipsoid is not cut away
  G4double pzBottomCut2 = 0 * um;  // here the bottom of the ellipsoid is cut
                                   // away below z = 0 (original configuration)
  // G4double pzBottomCut2 = 18.07*um; // test with position of the slice of
  // interest for vis.mac visualisation of the slice z = 18.07 µm

  // Ellipse2 = new G4Ellipsoid("Ellipse2",halfx2,halfy2,halfz2); // test with
  // entire ellipsoid
  ellipse2 = new G4Ellipsoid("Ellipse2", halfx2, halfy2, halfz2, pzBottomCut2,
                             pzTopCut2);  // cut ellipsoid

  logicEllipse2 = new G4LogicalVolume(ellipse2,  // its solid
                                      material2,  // its material
                                      "Ellipse2");  // its name

  physiEllipse2 = new G4PVPlacement(nullptr,  // no rotation
                                    center2,  // its position
                                    logicEllipse2,  // its logical volume
                                    "Ellipse2",  // its name
                                    logicEllipse1,  // its mother
                                    false,  // no boolean operation
                                    0);  // copy number

  // THIRD ELLIPSOID : Nucleus of cell 1 (contained inside ellipsoid 2)
  material3 = nist->FindOrBuildMaterial("Cell1_10times");
  G4ThreeVector center3 =
    G4ThreeVector(2.36 * um, -7.09 * um, 18.07 * um);  // position according to ellipsoid 2
  G4double halfx3 = 1.95 * um;  // Semiaxis in X
  G4double halfy3 = 3.23 * um;  // Semiaxis in Y
  G4double halfz3 = 4.32 * um;  // Semiaxis in Z

  // if we want to visualize slice at z = 18.07 µm using vis.mac,
  // we should cut all the ellipsoids (cut away all parts below 18.07 µm)
  // Visualization test: the lower part of ellipse 3 is cut off
  // G4double pzTopCut3 =4.32 *um;  // Visualization test
  // G4double pzBottomCut3 = 0*um;  // Visualization test: cut at zero before
  // moving ellipsoid 3
  // // to position defined just above => it will be cut at 18.07 µm after
  // moving

  ellipse3 = new G4Ellipsoid("Ellipse3", halfx3, halfy3, halfz3);  // entire ellipsoid
  // ellipse3 = new
  // G4Ellipsoid("Ellipse3",halfx3,halfy3,halfz3,pzBottomCut3,pzTopCut3); //
  // Visualization test

  logicEllipse3 = new G4LogicalVolume(ellipse3,  // its solid
                                      material3,  // its material
                                      "Ellipse3");  // its name

  physiEllipse3 = new G4PVPlacement(nullptr,  // no rotation
                                    center3,  // its position
                                    logicEllipse3,  // its logical volume
                                    "Ellipse3",  // its name
                                    logicEllipse2,  // its mother
                                    false,  // no boolean operation
                                    0);  // copy number

  // FOURTH ELLIPSOID : Nucleus of cell 2 (contained inside ellipsoid 2)
  material4 = nist->FindOrBuildMaterial("Cell2_10times");
  G4ThreeVector center4 =
    G4ThreeVector(8.66 * um, -3.15 * um, 18.07 * um);  // position according to ellipsoid 2
  G4double halfx4 = 2.08 * um;  // Semiaxis in X
  G4double halfy4 = 2.46 * um;  // Semiaxis in Y
  G4double halfz4 = 4.32 * um;  // Semiaxis in Z

  // means the lower part of ellipse 4 was cut off
  // G4double pzTopCut4 = 4.32*um;
  // G4double pzBottomCut4 = 0*um;// Visualization test: cut at zero before
  // moving ellipsoid 4
  // // to position defined just above => it will be cut at 18.07 µm after
  // moving

  ellipse4 = new G4Ellipsoid("Ellipse4", halfx4, halfy4, halfz4);  // entire ellipsoid
  // Ellipse4 = new
  // G4Ellipsoid("Ellipse4",halfx4,halfy4,halfz4,pzBottomCut4,pzTopCut4); //
  // Visualization test

  logicEllipse4 = new G4LogicalVolume(ellipse4,  // its solid
                                      material4,  // its material
                                      "Ellipse4");  // its name

  physiEllipse4 = new G4PVPlacement(nullptr,  // no rotation
                                    center4,  // its position
                                    logicEllipse4,  // its logical volume
                                    "Ellipse4",  // its name
                                    logicEllipse2,  // its mother
                                    false,  // no boolean operation
                                    0);  // copy number

  // FIFTH ELLIPSOID: intestine of the worm
  // Original (real) composition in dry mass (the worm is dehydrated by
  // lyophilization) Each element content is defined by its mass fraction (must
  // be normalised to 1)
  material5 = nist->FindOrBuildMaterial("Intestine_10times");
  G4ThreeVector center5 =
    G4ThreeVector(1.64 * um, 0.61 * um, 0 * um);  // position according to ellipsoid 2
  G4RotationMatrix rm5;
  G4double angle5 = -58.986 * deg;  // ellipse 5 is rotated around z axis
                                    // (rotation defined according to ellipse 2)
  // G4double angle5=0*deg;
  rm5.rotateZ(angle5);

  G4double halfx5 = 3.67 * um;  // Semiaxis in X
  G4double halfy5 = 16.48 * um;  // Semiaxis in Y
  G4double halfz5 = 28.68 * um;  // Semiaxis in Z

  G4double pzTopCut5 = 28.68 * um;  // here the top of the ellipsoid is not cut away
  G4double pzBottomCut5 = 0 * um;  // here the bottom of the ellipsoid is cut
                                   // away below z = 0 // original configuration
  // G4double pzBottomCut5 =18.07*um; // position of the slice of interest    //
  // Visualization test

  // ellipse5 = new G4Ellipsoid("Ellipse5",halfx5,halfy5,halfz5); // entire
  // ellipsoid // Visualization test
  ellipse5 = new G4Ellipsoid("Ellipse5", halfx5, halfy5, halfz5, pzBottomCut5,
                             pzTopCut5);  // original configuration

  logicEllipse5 = new G4LogicalVolume(ellipse5,  // its solid
                                      material5,  // its material
                                      "Ellipse5");  // its name

  physiEllipse5 = new G4PVPlacement(G4Transform3D(rm5, center5),  // rotation
                                    logicEllipse5,  // its logical volume
                                    "Ellipse5",  // its name
                                    logicEllipse2,  // its mother
                                    false,  // no boolean operation
                                    0);  // copy number

  // SIXTH ELLIPSOID: Nucleus of cell Ti (contained inside ellipsoid 5)
  // Original (real) composition in dry mass (the worm is dehydrated by
  // lyophilization) Each element content is defined by its mass fraction (must
  // be normalised to 1)

  material6 = nist->FindOrBuildMaterial("CellTi_1000times");
  G4ThreeVector center6 =
    G4ThreeVector(0.005 * um, 5.83 * um, 18.07 * um);  // position according to ellipsoid 5
  G4double halfx6 = 1.62 * um;  // Semiaxis in X
  G4double halfy6 = 1.95 * um;  // Semiaxis in Y
  G4double halfz6 = 1.62 * um;  // Semiaxis in Z

  // means the lower part of ellipse 6 was cut off
  // G4double pzTopCut6 = 1.62*um;
  // G4double pzBottomCut6 = 0*um;  // Visualization test: cut at zero before
  // moving ellipsoid 6
  // // to position defined just above => it will be cut at 18.07 µm after
  // moving. As ellipsoid 6 is contained in ellipsoid 5, it has been rotated at
  // the same time as ellipsoid 5. We want to counterbalance this rotation and
  // put ellipsoid 6 in the right direction (not rotated at all). So we must
  // rotate ellipsoid 6 of the opposite angle as ellipsoid 5
  G4RotationMatrix rm6;
  G4double angle6 = 58.986 * deg;  // ellipsoid 6 is rotated around z axis
                                   // (rotation defined according to ellipse 5)
  // G4double angle6= 0*deg;
  rm6.rotateZ(angle6);
  ellipse6 = new G4Ellipsoid("Ellipse6", halfx6, halfy6, halfz6);  // entire ellipsoid
  // ellipse6 = new
  // G4Ellipsoid("Ellipse6",halfx6,halfy6,halfz6,pzBottomCut6,pzTopCut6); 	 //
  // Visualization test
  logicEllipse6 = new G4LogicalVolume(ellipse6,  // its solid
                                      material6,  // its material
                                      "Ellipse6");  // its name

  physiEllipse6 =
    new G4PVPlacement(G4Transform3D(rm6, center6),  // So now ellipsoid 6 is vertical again
                                                    // (appears not rotated at all)
                      logicEllipse6,  // its logical volume
                      "Ellipse6",  // its name
                      logicEllipse5,  // its mother
                      false,  // no boolean operation
                      0);  // copy number

  /* logicEllipse1->SetUserLimits(new G4UserLimits(0.5*um,DBL_MAX,DBL_MAX,0));

  logicEllipse2->SetUserLimits(new G4UserLimits(0.5*um,DBL_MAX,DBL_MAX,0)); */

  // As the material composing the C. elegans worm is very transparent to
  // photons, the default steps in Geant4 (distance between pre-step and
  // post-step point) are too long. To perform a more precise simulation,
  // maximal values are defined for step length in each volume, according to its
  // size and configuration.

  logicEllipse1->SetUserLimits(new G4UserLimits(0.2 * um));
  logicEllipse2->SetUserLimits(new G4UserLimits(0.3 * um));
  logicEllipse3->SetUserLimits(new G4UserLimits(0.4 * um));
  logicEllipse4->SetUserLimits(new G4UserLimits(0.4 * um));
  logicEllipse5->SetUserLimits(new G4UserLimits(0.7 * um));
  logicEllipse6->SetUserLimits(new G4UserLimits(0.3 * um));
  // fLogicWorld->SetUserLimits(new G4UserLimits(1*um)); // world does not have
  // to be simulated as precisely

  //    fLogicWorld ->SetVisAttributes(new
  //    G4VisAttributes(G4Colour(1.0,1.0,1.0)));
  logicEllipse1->SetVisAttributes(new G4VisAttributes(G4Colour(1.0, 0.8, 1.0, 0.5)));
  logicEllipse2->SetVisAttributes(new G4VisAttributes(G4Colour(1.0, 0.5, 1.0, 0.5)));
  //    logicEllipse3 ->SetVisAttributes(new
  //    G4VisAttributes(G4Colour(0.0,1.0,1.0))); logicEllipse4
  //    ->SetVisAttributes(new G4VisAttributes(G4Colour(0.5,0.3,0.5)));
  logicEllipse5->SetVisAttributes(new G4VisAttributes(G4Colour(0.5, 0.5, 1.0, 0.5)));
  //    logicEllipse6 ->SetVisAttributes(new
  //    G4VisAttributes(G4Colour(0.5,0.8,0.2)));

  auto logVisAtt3 = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
  logVisAtt3->SetForceSolid(true);
  logicEllipse3->SetVisAttributes(logVisAtt3);

  auto logVisAtt4 = new G4VisAttributes(G4Colour(0.5, 0.3, 0.5));
  logVisAtt4->SetForceSolid(true);
  logicEllipse4->SetVisAttributes(logVisAtt4);

  auto logVisAtt6 = new G4VisAttributes(G4Colour(0.5, 0.8, 0.2));
  logVisAtt6->SetForceSolid(true);
  logicEllipse6->SetVisAttributes(logVisAtt6);

  PrintGeomParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::Construct_Phantom3()
{
  G4NistManager* nist = G4NistManager::Instance();
  material_GDP = nist->FindOrBuildMaterial("GDP");

  G4double pRmin = 171 * um;  // thickness = 25 um
  G4double pRmax = 196 * um;

  //    G4double pRmin = 146.8*um;  // thickness = 10 um
  //    G4double pRmax = 156.8*um;

  //    G4double pRmin = 293.6*um;  // thickness = 20 um
  //    G4double pRmax = 313.6*um;

  solid_GDP = new G4Sphere("GDP",  // name
                           pRmin, pRmax, 0., 2 * pi, 0., pi);

  logic_GDP = new G4LogicalVolume(solid_GDP,  // its solid
                                  material_GDP,  // its material
                                  "GDP");  // its name

  physi_GDP = new G4PVPlacement(nullptr,  // no rotation
                                G4ThreeVector(fXposAbs, 0., 0.),  // its position
                                logic_GDP,  // its logical volume
                                "GDP",  // its name
                                fLogicWorld,  // its mother
                                false,  // no boolean operation
                                0);  // copy number

  logic_GDP->SetUserLimits(new G4UserLimits(1 * um));

  auto logVisAtt = new G4VisAttributes(G4Colour(0.670588, 0.333333, 0.862745, 0.5));
  //    new G4VisAttributes(G4Colour(0.670588, 0.333333, 0.862745));
  logVisAtt->SetForceSolid(true);
  //  logVisAtt->SetVisibility(true);

  logic_GDP->SetVisAttributes(logVisAtt);

  PrintGeomParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

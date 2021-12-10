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
/// \file optical/wls/src/WLSDetectorConstruction.cc
/// \brief Implementation of the WLSDetectorConstruction class
//
//

#include "WLSDetectorConstruction.hh"

#include "WLSDetectorMessenger.hh"
#include "WLSMaterials.hh"
#include "WLSPhotonDetSD.hh"

#include "globals.hh"
#include "G4Box.hh"
#include "G4EllipticalTube.hh"
#include "G4ios.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4OpticalSurface.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSDetectorConstruction::WLSDetectorConstruction()
  : fVisAttributes()
  , fMaterials(nullptr)
  , fLogicHole(nullptr)
  , fLogicWorld(nullptr)
  , fPhysiWorld(nullptr)
  , fPhysiHole(nullptr)
{
  fDetectorMessenger = new WLSDetectorMessenger(this);

  fNumOfCladLayers = 0;

  fSurfaceRoughness = 1;

  fMirrorToggle       = true;
  fMirrorPolish       = 1.;
  fMirrorReflectivity = 1.;

  fMPPCPolish       = 1.;
  fMPPCReflectivity = 0.;

  fExtrusionPolish       = 1.;
  fExtrusionReflectivity = 1.;

  fXYRatio = 1.0;

  fWLSfiberZ      = 1. * m;
  fWLSfiberRY     = 0.5 * mm;
  fWLSfiberOrigin = 0.0;

  fMPPCShape = "Circle";
  fMPPCHalfL = fWLSfiberRY;
  fMPPCDist  = 0.00 * mm;
  fMPPCTheta = 0.0 * deg;
  fMPPCZ     = 0.05 * mm;

  fClrfiberZ = fMPPCZ + 10. * nm;
  fMirrorZ   = 0.1 * mm;

  fBarLength        = 1. * m;
  fBarBase          = 9.6 * mm;
  fHoleRadius       = 0.9 * mm;
  fHoleLength       = fBarLength;
  fCoatingThickness = 0.25 * mm;
  fCoatingRadius    = 1.875 * mm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSDetectorConstruction::~WLSDetectorConstruction()
{
  if(fDetectorMessenger)
    delete fDetectorMessenger;
  if(fMaterials)
    delete fMaterials;
  for (auto visAttributes: fVisAttributes)
  {
    delete visAttributes;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* WLSDetectorConstruction::Construct()
{
  if(fPhysiWorld)
  {
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();
    G4LogicalSkinSurface::CleanSurfaceTable();
    G4LogicalBorderSurface::CleanSurfaceTable();
  }

  fMaterials = WLSMaterials::GetInstance();
  UpdateGeometryParameters();

  return ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* WLSDetectorConstruction::ConstructDetector()
{

  auto air = FindMaterial("G4_AIR");
  //G4cout << "\nMaterial Properties Table for G4_AIR:" << G4endl;
  //air->GetMaterialPropertiesTable()->DumpTable();

  //--------------------------------------------------
  // World
  //--------------------------------------------------

  G4VSolid* solidWorld =
    new G4Box("World", fWorldSizeX, fWorldSizeY, fWorldSizeZ);

  fLogicWorld =
    new G4LogicalVolume(solidWorld, air, "World");

  fPhysiWorld =
    new G4PVPlacement(0, G4ThreeVector(), fLogicWorld, "World", 0, false, 0);

  //--------------------------------------------------
  // Extrusion
  //--------------------------------------------------

  auto coating = FindMaterial("Coating");

  G4VSolid* solidExtrusion = new G4Box("Extrusion", GetBarBase() / 2.,
                                       GetBarBase() / 2., GetBarLength() / 2.);

  G4LogicalVolume* logicExtrusion =
    new G4LogicalVolume(solidExtrusion, coating, "Extrusion");

  G4OpticalSurface* TiO2Surface = new G4OpticalSurface(
    "TiO2Surface", glisur, ground, dielectric_metal, fExtrusionPolish);

  G4MaterialPropertiesTable* TiO2SurfaceProperty =
    new G4MaterialPropertiesTable();

  std::vector<G4double> p_TiO2 = { 2.00 * eV, 3.47 * eV };

  std::vector<G4double> refl_TiO2 = { fExtrusionReflectivity,
                                      fExtrusionReflectivity };
  std::vector<G4double> effi_TiO2 = { 0., 0. };

  TiO2SurfaceProperty->AddProperty("REFLECTIVITY", p_TiO2, refl_TiO2);
  TiO2SurfaceProperty->AddProperty("EFFICIENCY", p_TiO2, effi_TiO2);

  TiO2Surface->SetMaterialPropertiesTable(TiO2SurfaceProperty);

  new G4PVPlacement(0, G4ThreeVector(), logicExtrusion, "Extrusion",
                    fLogicWorld, false, 0);

  new G4LogicalSkinSurface("TiO2Surface", logicExtrusion, TiO2Surface);

  //--------------------------------------------------
  // Scintillator
  //--------------------------------------------------

  auto polystyrene = FindMaterial("Polystyrene");
  //G4cout << "\nMaterial Properties Table for Polystyrene:" << G4endl;
  //polystyrene->GetMaterialPropertiesTable()->DumpTable();

  G4VSolid* solidScintillator =
    new G4Box("Scintillator",
              GetBarBase() / 2. - GetCoatingThickness() - GetCoatingRadius(),
              GetBarBase() / 2. - GetCoatingThickness() - GetCoatingRadius(),
              GetBarLength() / 2.);

  G4LogicalVolume* logicScintillator = new G4LogicalVolume(
    solidScintillator, polystyrene, "Scintillator");

  new G4PVPlacement(0, G4ThreeVector(), logicScintillator, "Scintillator",
                    logicExtrusion, false, 0);
  
  G4LogicalVolume* logicScintSide = nullptr;
  G4LogicalVolume* logicScintCrnr = nullptr;
  if(GetCoatingRadius() > 0.)
  {
    G4VSolid* solidScintside =
      new G4Box("SideOfBar",
                GetBarBase() / 2. - GetCoatingThickness() - GetCoatingRadius(),
                GetCoatingRadius() / 2., GetBarLength() / 2.);

    G4VSolid* solidScintcrnr =
      new G4Tubs("CrnrOfBar", 0.0 * cm, GetCoatingRadius(), GetBarLength() / 2.,
                 0. * deg, 90. * deg);

    logicScintSide = new G4LogicalVolume(
      solidScintside, polystyrene, "SideOfBar");

    logicScintCrnr = new G4LogicalVolume(
      solidScintcrnr, polystyrene, "CrnrOfBar");

    G4double pos =
      GetBarBase() / 2. - GetCoatingThickness() - GetCoatingRadius() / 2.;

    new G4PVPlacement(nullptr, G4ThreeVector(0., -pos, 0.), logicScintSide,
                      "SideOfBar", logicExtrusion, false, 0);

    new G4PVPlacement(nullptr, G4ThreeVector(0., pos, 0.), logicScintSide,
                      "SideOfBar", logicExtrusion, false, 1);

    G4RotationMatrix* rot1 = new G4RotationMatrix();
    rot1->rotateZ(-90.*deg);

    new G4PVPlacement(rot1, G4ThreeVector(pos, 0., 0.), logicScintSide,
                      "SideOfBar", logicExtrusion, false, 2);

    new G4PVPlacement(rot1, G4ThreeVector(-pos, 0., 0.), logicScintSide,
                      "SideOfBar", logicExtrusion, false, 3);

    pos = GetBarBase() / 2. - GetCoatingThickness() - GetCoatingRadius();

    new G4PVPlacement(nullptr, G4ThreeVector(pos, pos, 0.), logicScintCrnr,
                      "CrnrOfBar", logicExtrusion, false, 0);

    new G4PVPlacement(rot1, G4ThreeVector(-pos, pos, 0.), logicScintCrnr,
                      "CrnrOfBar", logicExtrusion, false, 1);

    G4RotationMatrix* rot2 = new G4RotationMatrix();
    rot2->rotateZ(-180.*deg);

    new G4PVPlacement(rot2, G4ThreeVector(-pos, -pos, 0.), logicScintCrnr,
                      "CrnrOfBar", logicExtrusion, false, 2);

    G4RotationMatrix* rot3 = new G4RotationMatrix();
    rot3->rotateZ(-270.*deg);

    new G4PVPlacement(rot3, G4ThreeVector(pos, -pos, 0.), logicScintCrnr,
                      "CrnrOfBar", logicExtrusion, false, 3);
  }

  if(GetFiberRadius() < GetHoleRadius())
  {
    G4VSolid* solidHole = new G4Tubs(
      "Hole", 0., GetHoleRadius(), GetHoleLength() / 2., 0. * deg, 360. * deg);

    fLogicHole = new G4LogicalVolume(solidHole, air, "Hole");

    fPhysiHole = new G4PVPlacement(0, G4ThreeVector(), fLogicHole, "Hole",
                                   logicScintillator, false, 0);
  }

  //--------------------------------------------------
  // Fiber
  //--------------------------------------------------

  if(!(fLogicHole) || !(fPhysiHole))
  {
    G4ExceptionDescription ed;
    ed << "The Fiber Hole has not been constructed";
    G4Exception("WLSDetectorConstruction", "wls001", FatalException, ed);
  }

  // Pointers to the most recently constructed volume
  G4LogicalVolume* logicPlacement   = fLogicHole;
  G4VPhysicalVolume* physiPlacement = fPhysiHole;

  //--------------------------------------------------
  // Fiber Construction
  //--------------------------------------------------

  // Boundary Surface Properties
  G4OpticalSurface* opSurface = nullptr;

  if(fSurfaceRoughness < 1.)
    opSurface = new G4OpticalSurface("RoughSurface", glisur, ground,
                                     dielectric_dielectric, fSurfaceRoughness);
  
  G4LogicalVolume* logicWLSfiber = nullptr;
  G4LogicalVolume* logicClad1    = nullptr;
  G4LogicalVolume* logicClad2    = nullptr;
  G4VPhysicalVolume* physiClad1  = nullptr;
  G4VPhysicalVolume* physiClad2  = nullptr;

  auto fpethylene = FindMaterial("FPethylene");
  auto pethylene = FindMaterial("Pethylene");
  auto pmma = FindMaterial("PMMA");

  // Determine the number of cladding layers to be built
  switch(fNumOfCladLayers)
  {
    case 2:

      //--------------------------------------------------
      // Cladding 2
      //--------------------------------------------------

      //G4cout << "\nMaterial Properties Table for fPethylene:" << G4endl;
      //fpethylene->GetMaterialPropertiesTable()->DumpTable();

      G4VSolid* solidClad2;

      if(fXYRatio == 1.)
        solidClad2 = new G4Tubs("Clad2", 0., fClad2RX, fClad2Z, 0., twopi);
      else
        solidClad2 = new G4EllipticalTube("Clad2", fClad2RX, fClad2RY, fClad2Z);

      logicClad2 =
        new G4LogicalVolume(solidClad2, fpethylene, "Clad2");

      physiClad2 =
        new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, fWLSfiberOrigin),
                          logicClad2, "Clad2", logicPlacement, false, 0);

      // Place the rough surface only if needed
      if(opSurface)
      {
        new G4LogicalBorderSurface("surfaceClad2Out", physiClad2,
                                   physiPlacement, opSurface);
        new G4LogicalBorderSurface("surfaceClad2In", physiPlacement, physiClad2,
                                   opSurface);
      }

      logicPlacement = logicClad2;
      physiPlacement = physiClad2;
      [[fallthrough]];

    case 1:

      //--------------------------------------------------
      // Cladding 1
      //--------------------------------------------------

      //G4cout << "\nMaterial Properties Table for Pethylene:" << G4endl;
      //pethylene->GetMaterialPropertiesTable()->DumpTable();

      G4VSolid* solidClad1;

      if(fXYRatio == 1.)
        solidClad1 = new G4Tubs("Clad1", 0., fClad1RX, fClad1Z, 0., twopi);
      else
        solidClad1 = new G4EllipticalTube("Clad1", fClad1RX, fClad1RY, fClad1Z);

      logicClad1 =
        new G4LogicalVolume(solidClad1, pethylene, "Clad1");

      physiClad1 =
        new G4PVPlacement(0, G4ThreeVector(0., 0., fWLSfiberOrigin), logicClad1,
                          "Clad1", logicPlacement, false, 0);

      // Place the rough surface only if needed
      if(opSurface)
      {
        new G4LogicalBorderSurface("surfaceClad1Out", physiClad1,
                                   physiPlacement, opSurface);

        new G4LogicalBorderSurface("surfaceClad1In", physiPlacement, physiClad1,
                                   opSurface);
      }

      logicPlacement = logicClad1;
      physiPlacement = physiClad1;
      [[fallthrough]];

    default:

      //--------------------------------------------------
      // WLS Fiber
      //--------------------------------------------------

      //G4cout << "\nMaterial Properties Table for PMMA:" << G4endl;
      //pmma->GetMaterialPropertiesTable()->DumpTable();

      G4VSolid* solidWLSfiber;

      if(fXYRatio == 1.)
      {
        solidWLSfiber =
          new G4Tubs("WLSFiber", 0., fWLSfiberRX, fWLSfiberZ, 0., twopi);
      }
      else
      {
        solidWLSfiber = new G4EllipticalTube("WLSFiber", fWLSfiberRX,
                                             fWLSfiberRY, fWLSfiberZ);
      }

      logicWLSfiber =
        new G4LogicalVolume(solidWLSfiber, pmma, "WLSFiber");

      logicWLSfiber->SetUserLimits(
        new G4UserLimits(DBL_MAX, DBL_MAX, 10. * ms));

      G4VPhysicalVolume* physiWLSfiber =
        new G4PVPlacement(0, G4ThreeVector(0., 0., fWLSfiberOrigin),
                          logicWLSfiber, "WLSFiber", logicPlacement, false, 0);

      // Place the rough surface only if needed
      if(opSurface)
      {
        new G4LogicalBorderSurface("surfaceWLSOut", physiWLSfiber,
                                   physiPlacement, opSurface);

        new G4LogicalBorderSurface("surfaceWLSIn", physiPlacement,
                                   physiWLSfiber, opSurface);
      }
  }

  //--------------------------------------------------
  // Mirror for reflection at one of the end
  //--------------------------------------------------

  // Place the mirror only if the user wants the mirror
  G4LogicalVolume* logicMirror = nullptr;

  auto aluminum = FindMaterial("G4_Al");

  if(fMirrorToggle)
  {
    G4VSolid* solidMirror =
      new G4Box("Mirror", fMirrorRmax, fMirrorRmax, fMirrorZ);

    logicMirror =
      new G4LogicalVolume(solidMirror, aluminum, "Mirror");

    G4OpticalSurface* mirrorSurface = new G4OpticalSurface(
      "MirrorSurface", glisur, ground, dielectric_metal, fMirrorPolish);

    G4MaterialPropertiesTable* mirrorSurfaceProperty =
      new G4MaterialPropertiesTable();

    std::vector<G4double> p_mirror    = { 2.00 * eV, 3.47 * eV };
    std::vector<G4double> refl_mirror = { fMirrorReflectivity,
                                          fMirrorReflectivity };
    std::vector<G4double> effi_mirror = { 0., 0. };

    mirrorSurfaceProperty->AddProperty("REFLECTIVITY", p_mirror, refl_mirror);
    mirrorSurfaceProperty->AddProperty("EFFICIENCY", p_mirror, effi_mirror);

    mirrorSurface->SetMaterialPropertiesTable(mirrorSurfaceProperty);

    new G4PVPlacement(0, G4ThreeVector(0., 0., fMirrorOrigin), logicMirror,
                      "Mirror", fLogicWorld, false, 0);

    new G4LogicalSkinSurface("MirrorSurface", logicMirror, mirrorSurface);
  }

  //--------------------------------------------------
  // Coupling at the read-out end
  //--------------------------------------------------

  // Clear Fiber (Coupling Layer)
  G4VSolid* solidCouple = new G4Box("Couple", fCoupleRX, fCoupleRY, fCoupleZ);

  G4LogicalVolume* logicCouple =
    new G4LogicalVolume(solidCouple, air, "Couple");

  new G4PVPlacement(0, G4ThreeVector(0., 0., fCoupleOrigin), logicCouple,
                    "Couple", fLogicWorld, false, 0);

  //--------------------------------------------------
  // A logical layer in front of PhotonDet
  //--------------------------------------------------

  // Purpose: Preventing direct dielectric to metal contact

  // Check for valid placement of PhotonDet
  if(fMPPCTheta > std::atan(fMPPCDist / fMPPCHalfL))
  {
    fMPPCTheta   = 0.;
    fMPPCOriginX = std::sin(fMPPCTheta) * (fMPPCDist + fClrfiberZ);
    fMPPCOriginZ = -fCoupleZ + std::cos(fMPPCTheta) * (fMPPCDist + fClrfiberZ);
    G4ExceptionDescription ed;
    ed << "Invalid alignment.  Alignment reset to 0.";
    G4Exception("WLSDetectorConstruction", "wls002", JustWarning, ed);
  }

  // Clear Fiber (Coupling Layer)
  G4VSolid* solidClrfiber;

  if(fMPPCShape == "Square")
  {
    solidClrfiber =
      new G4Box("ClearFiber", fClrfiberHalfL, fClrfiberHalfL, fClrfiberZ);
  }
  else
  {
    solidClrfiber =
      new G4Tubs("ClearFiber", 0., fClrfiberHalfL, fClrfiberZ, 0., twopi);
  }

  G4LogicalVolume* logicClrfiber =
    new G4LogicalVolume(solidClrfiber, air, "ClearFiber");

  new G4PVPlacement(new G4RotationMatrix(CLHEP::HepRotationY(-fMPPCTheta)),
                    G4ThreeVector(fMPPCOriginX, 0.0, fMPPCOriginZ),
                    logicClrfiber, "ClearFiber", logicCouple, false, 0);

  //--------------------------------------------------
  // PhotonDet (Sensitive Detector)
  //--------------------------------------------------

  // Physical Construction
  G4VSolid* solidPhotonDet = nullptr;

  if(fMPPCShape == "Square")
    solidPhotonDet = new G4Box("PhotonDet", fMPPCHalfL, fMPPCHalfL, fMPPCZ);
  else
    solidPhotonDet = new G4Tubs("PhotonDet", 0., fMPPCHalfL, fMPPCZ, 0., twopi);

  G4LogicalVolume* logicPhotonDet =
    new G4LogicalVolume(solidPhotonDet, aluminum, "PhotonDet_LV");

  new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicPhotonDet, "PhotonDet",
                    logicClrfiber, false, 0);

  // PhotonDet Surface Properties
  G4OpticalSurface* photonDetSurface = new G4OpticalSurface(
    "PhotonDetSurface", glisur, ground, dielectric_metal, fMPPCPolish);

  G4MaterialPropertiesTable* photonDetSurfaceProperty =
    new G4MaterialPropertiesTable();

  std::vector<G4double> p_mppc    = { 2.00 * eV, 3.47 * eV };
  std::vector<G4double> refl_mppc = { fMPPCReflectivity, fMPPCReflectivity };
  std::vector<G4double> effi_mppc = { 1., 1. };

  photonDetSurfaceProperty->AddProperty("REFLECTIVITY", p_mppc, refl_mppc);
  photonDetSurfaceProperty->AddProperty("EFFICIENCY", p_mppc, effi_mppc);

  photonDetSurface->SetMaterialPropertiesTable(photonDetSurfaceProperty);

  new G4LogicalSkinSurface("PhotonDetSurface", logicPhotonDet,
                           photonDetSurface);

  // visualization attributes -------------------------------------------------

  auto visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  visAttributes->SetVisibility(false);
  fLogicWorld->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour(0.2,0.2,0.2,0.5));
  visAttributes->SetVisibility(true);
  logicExtrusion->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour(0.0,0.0,1.0,0.9));
  visAttributes->SetVisibility(true);
  logicScintillator->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour(0.0,0.8,0.2,0.2));
  visAttributes->SetVisibility(true);
  logicScintSide->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour(0.0,0.8,0.2,0.2));
  visAttributes->SetVisibility(true);
  logicScintCrnr->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour(0.4,0.0,0.0,0.5));
  visAttributes->SetVisibility(true);
  fLogicHole->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  if(logicClad1 != nullptr)
  {
    visAttributes = new G4VisAttributes(G4Colour(0.0,0.8,0.5,0.5));
    visAttributes->SetVisibility(true);
    logicClad1->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
  }

  if(logicClad2 != nullptr)
  {
    visAttributes = new G4VisAttributes(G4Colour(0.0,0.5,0.8,0.5));
    visAttributes->SetVisibility(true);
    logicClad2->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
  }

  visAttributes = new G4VisAttributes(G4Colour(0.8,0.8,1.0));
  visAttributes->SetVisibility(true);
  logicWLSfiber->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  if(fMirrorToggle == true)
  {
    visAttributes = new G4VisAttributes(G4Colour(0.3,0.3,1.0,0.3));
    visAttributes->SetVisibility(true);
    logicMirror->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
  }

  visAttributes = new G4VisAttributes(G4Colour(0.0,0.0,0.5,0.5));
  visAttributes->SetVisibility(true);
  logicCouple->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour(0.3,0.3,0.3,0.5));
  visAttributes->SetVisibility(true);
  logicClrfiber->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0,0.8));
  visAttributes->SetVisibility(true);
  logicPhotonDet->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  return fPhysiWorld;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::ConstructSDandField()
{
  if(!fmppcSD.Get())
  {
    G4String mppcSDName    = "WLS/PhotonDet";
    WLSPhotonDetSD* mppcSD = new WLSPhotonDetSD(mppcSDName);
    G4SDManager::GetSDMpointer()->AddNewDetector(mppcSD);
    fmppcSD.Put(mppcSD);
  }
  SetSensitiveDetector("PhotonDet_LV", fmppcSD.Get(), true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::UpdateGeometryParameters()
{
  fWLSfiberRX = fXYRatio * fWLSfiberRY;

  fClad1RX = fWLSfiberRX + 0.03 * fWLSfiberRX;
  fClad1RY = fWLSfiberRY + 0.03 * fWLSfiberRY;
  fClad1Z  = fWLSfiberZ;

  fClad2RX = fClad1RX + 0.03 * fWLSfiberRX;
  fClad2RY = fClad1RY + 0.03 * fWLSfiberRY;
  fClad2Z  = fWLSfiberZ;

  fWorldSizeX = fClad2RX + fMPPCDist + fMPPCHalfL + 1. * cm;
  fWorldSizeY = fClad2RY + fMPPCDist + fMPPCHalfL + 1. * cm;
  fWorldSizeZ = fWLSfiberZ + fMPPCDist + fMPPCHalfL + 1. * cm;

  fCoupleRX = fWorldSizeX;
  fCoupleRY = fWorldSizeY;
  fCoupleZ  = (fWorldSizeZ - fWLSfiberZ) / 2.;

  fClrfiberHalfL = fMPPCHalfL;

  fMirrorRmax = fClad2RY;

  fCoupleOrigin = fWLSfiberOrigin + fWLSfiberZ + fCoupleZ;
  fMirrorOrigin = fWLSfiberOrigin - fWLSfiberZ - fMirrorZ;
  fMPPCOriginX  = std::sin(fMPPCTheta) * (fMPPCDist + fClrfiberZ);
  fMPPCOriginZ  = -fCoupleZ + std::cos(fMPPCTheta) * (fMPPCDist + fClrfiberZ);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetPhotonDetGeometry(G4String shape)
// Set the Geometry of the PhotonDet detector
// Pre:  shape must be either "Circle" and "Square"
{
  if(shape == "Circle" || shape == "Square")
    fMPPCShape = shape;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetNumberOfCladding(G4int num)
// Set the number of claddings
// Pre: 0 <= num <= 2
{
  fNumOfCladLayers = num;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetWLSLength(G4double length)
// Set the TOTAL length of the WLS fiber
{
  fWLSfiberZ = length;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetWLSRadius(G4double radius)
// Set the Y radius of WLS fiber
{
  fWLSfiberRY = radius;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetClad1Radius(G4double radius)
// Set the Y radius of Cladding 1
{
  fClad1RY = radius;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetClad2Radius(G4double radius)
// Set the Y radius of Cladding 2
{
  fClad2RY = radius;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetPhotonDetHalfLength(G4double halfL)
// Set the half length of the PhotonDet detector
// The half length will be the radius if PhotonDet is circular
{
  fMPPCHalfL = halfL;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetGap(G4double gap)
// Set the distance between fiber end and PhotonDet
{
  fMPPCDist = gap;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetPhotonDetAlignment(G4double theta)
// Set the Aligment of PhotonDet with respect to the z axis
// If theta is 0 deg, then the detector is perfectly aligned
// PhotonDet will be deviated by theta from z axis
// facing towards the center of the fiber
{
  fMPPCTheta = theta;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetSurfaceRoughness(G4double roughness)
// Set the Surface Roughness between Cladding 1 and WLS fiber
// Pre: 0 < roughness <= 1
{
  fSurfaceRoughness = roughness;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetMirrorPolish(G4double polish)
// Set the Polish of the mirror, polish of 1 is a perfect mirror surface
// Pre: 0 < polish <= 1
{
  fMirrorPolish = polish;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetMirrorReflectivity(G4double reflectivity)
// Set the Reflectivity of the mirror, reflectivity of 1 is a perfect mirror
// Pre: 0 < reflectivity <= 1
{
  fMirrorReflectivity = reflectivity;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetPhotonDetPolish(G4double polish)
// Set the Polish of the PhotonDet, polish of 1 is a perfect mirror surface
// Pre: 0 < polish <= 1
{
  fMPPCPolish = polish;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetPhotonDetReflectivity(G4double reflectivity)
// Set the Reflectivity of the PhotonDet, reflectivity of 1 is a perfect mirror
// Pre: 0 < reflectivity <= 1
{
  fMPPCReflectivity = reflectivity;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetMirror(G4bool flag)
// Toggle to place the mirror or not at one end (-z end) of the fiber
// True means place the mirror, false means otherwise
{
  fMirrorToggle = flag;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetXYRatio(G4double r)
// Set the ratio of the x and y radius of the ellipse (x/y)
// a ratio of 1 would produce a circle
{
  fXYRatio = r;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetBarLength(G4double length)
// Set the length of the scintillator bar
{
  fBarLength = length;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetBarBase(G4double side)
// Set the side of the scintillator bar
{
  fBarBase = side;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetHoleRadius(G4double radius)
// Set the radius of the fiber hole
{
  fHoleRadius = radius;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetCoatingThickness(G4double thick)
// Set thickness of the coating on the bars
{
  fCoatingThickness = thick;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetCoatingRadius(G4double radius)
// Set inner radius of the corner bar coating
{
  fCoatingRadius = radius;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double WLSDetectorConstruction::GetWLSFiberLength() { return fWLSfiberZ; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double WLSDetectorConstruction::GetBarLength() { return fBarLength; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double WLSDetectorConstruction::GetBarBase() { return fBarBase; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double WLSDetectorConstruction::GetHoleRadius() { return fHoleRadius; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double WLSDetectorConstruction::GetHoleLength() { return fHoleLength; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double WLSDetectorConstruction::GetFiberRadius() { return GetWLSFiberRMax(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double WLSDetectorConstruction::GetCoatingThickness()
{
  return fCoatingThickness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double WLSDetectorConstruction::GetCoatingRadius() { return fCoatingRadius; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double WLSDetectorConstruction::GetWLSFiberEnd()
{
  return fWLSfiberOrigin + fWLSfiberZ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double WLSDetectorConstruction::GetWLSFiberRMax()
{
  if(fNumOfCladLayers == 2)
    return fClad2RY;
  if(fNumOfCladLayers == 1)
    return fClad1RY;
  return fWLSfiberRY;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double WLSDetectorConstruction::GetSurfaceRoughness()
{
  return fSurfaceRoughness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Return True if the fiber construction is ideal
G4bool WLSDetectorConstruction::IsPerfectFiber()
{
  return fSurfaceRoughness == 1. && fXYRatio == 1. &&
         (!fMirrorToggle || (fMirrorPolish == 1. && fMirrorReflectivity == 1.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* WLSDetectorConstruction::FindMaterial(G4String name)
{
  G4Material* material = G4Material::GetMaterial(name, true);
  return material;
}

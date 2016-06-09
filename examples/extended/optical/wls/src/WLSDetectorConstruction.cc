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

#include "G4ios.hh"
#include "globals.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4OpBoundaryProcess.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4GeometryManager.hh"
#include "G4SDManager.hh"

#include "G4SolidStore.hh"
#include "G4RegionStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4RunManager.hh"

#include "WLSDetectorConstruction.hh"
#include "WLSDetectorMessenger.hh"
#include "WLSMaterials.hh"
#include "WLSPhotonDetSD.hh"

#include "G4UserLimits.hh"

WLSPhotonDetSD* WLSDetectorConstruction::mppcSD = NULL;

WLSDetectorConstruction::WLSDetectorConstruction()
 : physiWorld(NULL)
{

  detectorMessenger = new WLSDetectorMessenger(this);
  materials = NULL;

  numOfCladLayers = 0;
 
  surfaceRoughness = 1;
 
  mirrorToggle = true;
  mirrorPolish = 1.;
  mirrorReflectivity = 1.;

  mppcPolish = 1.;
  mppcReflectivity = 0.;

  extrusionPolish = 1.;
  extrusionReflectivity = 1.;
 
  XYRatio = 1.0;

  wlsfiberZ     = 1.*m;
  wlsfiberRY  = 0.5*mm;
  wlsfiberOrigin = 0.0;
 
  mppcShape = "Circle";
  mppcHalfL = wlsfiberRY;
  mppcDist  = 0.00*mm;
  mppcTheta = 0.0*deg;
  mppcZ     = 0.05*mm;
 
  clrfiberZ  = mppcZ + 10.*nm ;
  mirrorZ    = 0.1*mm;

  barLength        = 1.*m;
  barBase          = 9.6*mm;
  holeRadius       = 0.9*mm;
  holeLength       = barLength;
  coatingThickness = 0.25*mm;
  coatingRadius    = 1.875*mm;

  UpdateGeometryParameters();
}

WLSDetectorConstruction::~WLSDetectorConstruction()
{
  if (detectorMessenger) delete detectorMessenger;
  if (materials)         delete materials;
}

G4VPhysicalVolume* WLSDetectorConstruction::Construct()
{
  materials = WLSMaterials::GetInstance();

  return ConstructDetector();
}

G4VPhysicalVolume* WLSDetectorConstruction::ConstructDetector()
{
  //--------------------------------------------------
  // World
  //--------------------------------------------------

  G4VSolid* solidWorld =
                       new G4Box("World", worldSizeX, worldSizeY, worldSizeZ);

  logicWorld = new G4LogicalVolume(solidWorld,
                                   FindMaterial("G4_AIR"),
                                   "World");

  physiWorld = new G4PVPlacement(0,
                                 G4ThreeVector(),
                                 logicWorld,
                                 "World",
                                 0,
                                 false,
                                 0);

  //--------------------------------------------------
  // Extrusion
  //--------------------------------------------------

  G4VSolid* solidExtrusion =
        new G4Box("Extrusion",GetBarBase()/2,GetBarBase()/2,GetBarLength()/2);

  G4LogicalVolume* logicExtrusion =
                      new G4LogicalVolume(solidExtrusion,
                                          FindMaterial("Coating"),
                                          "Extrusion");

  G4OpticalSurface* TiO2Surface = new G4OpticalSurface("TiO2Surface",
                                                       glisur,
                                                       ground,
                                                       dielectric_metal,
                                                       extrusionPolish);

  G4MaterialPropertiesTable* TiO2SurfaceProperty =
                                             new G4MaterialPropertiesTable();

  G4double p_TiO2[2] = {2.00*eV, 3.47*eV};
  G4double refl_TiO2[2] = {extrusionReflectivity,extrusionReflectivity};
  G4double effi_TiO2[2] = {0, 0};

  TiO2SurfaceProperty -> AddProperty("REFLECTIVITY",p_TiO2,refl_TiO2,2);
  TiO2SurfaceProperty -> AddProperty("EFFICIENCY",p_TiO2,effi_TiO2,2);

  TiO2Surface -> SetMaterialPropertiesTable(TiO2SurfaceProperty);

  new G4PVPlacement(0,
                    G4ThreeVector(),
                    logicExtrusion,
                    "Extrusion",
                    logicWorld,
                    false,
                    0);

  new G4LogicalSkinSurface("TiO2Surface",logicExtrusion,TiO2Surface);

  //--------------------------------------------------
  // Scintillator
  //--------------------------------------------------

  G4VSolid* solidScintillator = new G4Box("Scintillator",
                                GetBarBase()/2-GetCoatingThickness()
                                           -GetCoatingRadius(),
                                GetBarBase()/2-GetCoatingThickness()
                                           -GetCoatingRadius(),
                                GetBarLength()/2);

  G4LogicalVolume* logicScintillator =
                             new G4LogicalVolume(solidScintillator,
                                                 FindMaterial("Polystyrene"),
                                                 "Scintillator");

  new G4PVPlacement(0,
                    G4ThreeVector(),
                    logicScintillator,
                    "Scintillator",
                    logicExtrusion,
                    false,
                    0);

  if (GetCoatingRadius() > 0.*mm) {
     G4VSolid* solidScintside = new G4Box("SideOfBar",
                                GetBarBase()/2-GetCoatingThickness()
                                           -GetCoatingRadius(),
                                GetCoatingRadius()/2,
                                GetBarLength()/2);
     G4VSolid* solidScintcrnr = new G4Tubs("CrnrOfBar",
                                 0.0*cm,
                                 GetCoatingRadius(),
                                 GetBarLength()/2,
                                 0.*deg,
                                 90.*deg);
     G4LogicalVolume* logicScintSide =
                             new G4LogicalVolume(solidScintside,
                                                 FindMaterial("Polystyrene"),
                                                 "SideOfBar");

     G4LogicalVolume* logicScintCrnr =
                             new G4LogicalVolume(solidScintcrnr,
                                                 FindMaterial("Polystyrene"),
                                                 "CrnrOfBar");

     G4double x=GetBarBase()/2-GetCoatingThickness()-GetCoatingRadius()/2;
     G4double y=GetBarBase()/2-GetCoatingThickness()-GetCoatingRadius()/2;

     new G4PVPlacement(0,
                       G4ThreeVector(0,-y,0),
                       logicScintSide,
                       "SideOfBar",
                       logicExtrusion,
                       false,
                       0);
     new G4PVPlacement(0,
                       G4ThreeVector(0, y,0),
                       logicScintSide,
                       "SideOfBar",
                       logicExtrusion,
                       false,
                       1);

     G4RotationMatrix* g4rot = new G4RotationMatrix();
     *g4rot = stringToRotationMatrix("Z90");
     *g4rot = g4rot->inverse();
     if (*g4rot == G4RotationMatrix()) g4rot = NULL;

     new G4PVPlacement(g4rot,
                       G4ThreeVector(x,0,0),
                       logicScintSide,
                       "SideOfBar",
                       logicExtrusion,
                       false,
                       2);
     new G4PVPlacement(g4rot,
                       G4ThreeVector(-x,0,0),
                       logicScintSide,
                       "SideOfBar",
                       logicExtrusion,
                       false,
                       3);

     x = GetBarBase()/2-GetCoatingThickness()-GetCoatingRadius();
     y = GetBarBase()/2-GetCoatingThickness()-GetCoatingRadius();

     new G4PVPlacement(0,
                       G4ThreeVector(x,y,0),
                       logicScintCrnr,
                       "CrnrOfBar",
                       logicExtrusion,
                       false,
                       0);

     new G4PVPlacement(g4rot,
                       G4ThreeVector(-x,y,0),
                       logicScintCrnr,
                       "CrnrOfBar",
                       logicExtrusion,
                       false,
                       1);

     g4rot = new G4RotationMatrix();
     *g4rot = stringToRotationMatrix("Z180");
     *g4rot = g4rot->inverse();
     if (*g4rot == G4RotationMatrix()) g4rot = NULL;

     new G4PVPlacement(g4rot,
                       G4ThreeVector(-x,-y,0),
                       logicScintCrnr,
                       "CrnrOfBar",
                       logicExtrusion,
                       false,
                       2);

     g4rot = new G4RotationMatrix();
     *g4rot = stringToRotationMatrix("Z270");
     *g4rot = g4rot->inverse();
     if (*g4rot == G4RotationMatrix()) g4rot = NULL;

     new G4PVPlacement(g4rot,
                       G4ThreeVector(x,-y,0),
                       logicScintCrnr,
                       "CrnrOfBar",
                       logicExtrusion,
                       false,
                       3); 

  }

  if (GetFiberRadius()<GetHoleRadius()) {

        G4VSolid* solidHole = new G4Tubs("Hole",
                               0.0*cm,
                               GetHoleRadius(),
                               GetHoleLength()/2,
                               0.*deg,
                               360.*deg);
        logicHole = new G4LogicalVolume(solidHole,
                                        FindMaterial("G4_AIR"),
                                        "Hole");

        physiHole = new G4PVPlacement(0,
                                      G4ThreeVector(),
                                      logicHole,
                                      "Hole",
                                      logicScintillator,
                                      false,
                                      0);
  }

  //--------------------------------------------------
  // Fiber
  //--------------------------------------------------

  ConstructFiber();

  //--------------------------------------------------
  // End of Construction
  //--------------------------------------------------

  return physiWorld;
}

void WLSDetectorConstruction::ConstructFiber()
{
  if (!(logicHole) || !(physiHole) ) {
     std::ostringstream o;
     o << "The Fiber Hole has not been constructed";
     G4Exception("WLSDetectorConstruction::ConstructFiber","",
                  FatalException,o.str().c_str());
  }

  // Pointers to the most recently constructed volume
  G4LogicalVolume* logicPlacement = logicHole;
  G4VPhysicalVolume* physiPlacement = physiHole;

  //--------------------------------------------------
  // Fiber Construction
  //-------------------------------------------------- 

  // Boundary Surface Properties
  G4OpticalSurface* OpSurface = NULL;
 
  if (surfaceRoughness < 1.)
     OpSurface = new G4OpticalSurface("RoughSurface",          // Surface Name
                                      glisur,                  // SetModel
                                      ground,                  // SetFinish
                                      dielectric_dielectric,   // SetType
                                      surfaceRoughness);       // SetPolish

  G4LogicalVolume   *logicClad1, *logicClad2;
  G4VPhysicalVolume *physiClad1, *physiClad2;

  // Determine the number of cladding layers to be built
  switch ( numOfCladLayers ) {
 
    case 2:

     //--------------------------------------------------
     // Cladding 2
     //--------------------------------------------------

     G4VSolid* solidClad2;
 
     if (XYRatio == 1.)
       solidClad2 = new G4Tubs("Clad2",0.,clad2RX,clad2Z,0.0*rad,twopi*rad);
     else
       solidClad2 = new G4EllipticalTube("Clad2",clad2RX,clad2RY,clad2Z);

     logicClad2  = new G4LogicalVolume(solidClad2,
                                       FindMaterial("FPethylene"),
                                       "Clad2");

     physiClad2 = new G4PVPlacement(0,
                                    G4ThreeVector(0.0,0.0,wlsfiberOrigin),
                                    logicClad2,
                                    "Clad2",
                                    logicPlacement,
                                    false,
                                    0);

     // Place the rough surface only if needed
     if (OpSurface) {
       new G4LogicalBorderSurface("surfaceClad2Out",
                                  physiClad2,
                                  physiPlacement,
                                  OpSurface);
       new G4LogicalBorderSurface("surfaceClad2In",
                                  physiPlacement,
                                  physiClad2,
                                  OpSurface);
     }

     logicPlacement = logicClad2;
     physiPlacement = physiClad2;

    case 1:

     //--------------------------------------------------
     // Cladding 1
     //--------------------------------------------------

     G4VSolid* solidClad1;

     if (XYRatio == 1.)
       solidClad1 = new G4Tubs("Clad1",0.,clad1RX,clad1Z,0.0*rad,twopi*rad);
     else
       solidClad1 = new G4EllipticalTube("Clad1",clad1RX,clad1RY,clad1Z);

     logicClad1 = new G4LogicalVolume(solidClad1,
                                      FindMaterial("Pethylene"),
                                      "Clad1");

     physiClad1 = new G4PVPlacement(0,
                                    G4ThreeVector(0.0,0.0,wlsfiberOrigin),
                                    logicClad1,
                                    "Clad1",
                                    logicPlacement,
                                    false,
                                    0);

     // Place the rough surface only if needed
     if (OpSurface) {
       new G4LogicalBorderSurface("surfaceClad1Out",
                                  physiClad1,
                                  physiPlacement,
                                  OpSurface);
       new G4LogicalBorderSurface("surfaceClad1In",
                                  physiPlacement,
                                  physiClad1,
                                  OpSurface);
     }

     logicPlacement = logicClad1;
     physiPlacement = physiClad1;

    default:

     //--------------------------------------------------
     // WLS Fiber
     //--------------------------------------------------

     G4VSolid* solidWLSfiber;

     if (XYRatio == 1.)
       solidWLSfiber =
             new G4Tubs("WLSFiber",0.,wlsfiberRX,wlsfiberZ,0.0*rad,twopi*rad);
     else
       solidWLSfiber =
             new G4EllipticalTube("WLSFiber",wlsfiberRX,wlsfiberRY,wlsfiberZ);

     G4LogicalVolume*   logicWLSfiber =
                                     new G4LogicalVolume(solidWLSfiber,
                                                         FindMaterial("PMMA"),
                                                         "WLSFiber");

     logicWLSfiber->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,10*ms));

     G4VPhysicalVolume* physiWLSfiber = new G4PVPlacement(0,
                                       G4ThreeVector(0.0,0.0,wlsfiberOrigin),
                                       logicWLSfiber,
                                       "WLSFiber",
                                       logicPlacement,
                                       false,
                                       0);

     // Place the rough surface only if needed
     if (OpSurface) {
       new G4LogicalBorderSurface("surfaceWLSOut",
                                  physiWLSfiber,
                                  physiPlacement,
                                  OpSurface);
       new G4LogicalBorderSurface("surfaceWLSIn",
                                  physiPlacement,
                                  physiWLSfiber,
                                  OpSurface);
     }
  }

  //--------------------------------------------------
  // Mirror for reflection at one of the end
  //--------------------------------------------------

  // Place the mirror only if the user wants the mirror
  if (mirrorToggle) {  

     G4VSolid* solidMirror = new G4Box("Mirror",
                                       mirrorRmax,
                                       mirrorRmax,
                                       mirrorZ);
 
     G4LogicalVolume* logicMirror = new G4LogicalVolume(solidMirror,
                                                        FindMaterial("G4_Al"),
                                                        "Mirror");

     G4OpticalSurface* MirrorSurface = new G4OpticalSurface("MirrorSurface",
                                                             glisur,
                                                             ground,
                                                             dielectric_metal,
                                                             mirrorPolish);

     G4MaterialPropertiesTable* MirrorSurfaceProperty =
                                              new G4MaterialPropertiesTable();

     G4double p_mirror[2] = {2.00*eV, 3.47*eV};
     G4double refl_mirror[2] = {mirrorReflectivity,mirrorReflectivity};
     G4double effi_mirror[2] = {0, 0};

     MirrorSurfaceProperty->AddProperty("REFLECTIVITY",p_mirror,refl_mirror,2);
     MirrorSurfaceProperty->AddProperty("EFFICIENCY",p_mirror,effi_mirror,2);

     MirrorSurface -> SetMaterialPropertiesTable(MirrorSurfaceProperty);

     new G4PVPlacement(0,
                       G4ThreeVector(0.0,0.0,mirrorOrigin),
                       logicMirror,
                       "Mirror",
                       logicWorld,
                       false,
                       0);

     new G4LogicalSkinSurface("MirrorSurface",logicMirror,MirrorSurface);
  }

  //--------------------------------------------------
  // Coupling at the read-out end
  //--------------------------------------------------  

  // Clear Fiber (Coupling Layer)
  G4VSolid* solidCouple = new G4Box("Couple",coupleRX,coupleRY,coupleZ);

  G4LogicalVolume*   logicCouple = new G4LogicalVolume(solidCouple,
                                                       FindMaterial("G4_AIR"),
                                                       "Couple");

  new G4PVPlacement(0,
                    G4ThreeVector(0.0,0.0,coupleOrigin),
                    logicCouple,
                    "Couple",
                    logicWorld,
                    false,
                    0);

  //--------------------------------------------------
  // A logical layer in front of PhotonDet
  //--------------------------------------------------  

  // Purpose: Preventing direct dielectric to metal contact  

  // Check for valid placement of PhotonDet
  if (mppcTheta > std::atan(mppcDist / mppcHalfL)) {

     mppcTheta = 0;
     mppcOriginX  = std::sin(mppcTheta) * (mppcDist + clrfiberZ);
     mppcOriginZ  = -coupleZ + std::cos(mppcTheta) * (mppcDist + clrfiberZ);
     G4cerr << "Invalid alignment.  Alignment Reset to 0" << G4endl;     
  }
 
  // Clear Fiber (Coupling Layer)
  G4VSolid* solidClrfiber;
 
  if ( mppcShape == "Square" )
    solidClrfiber = 
                new G4Box("ClearFiber",clrfiberHalfL,clrfiberHalfL,clrfiberZ);
  else
    solidClrfiber =
       new G4Tubs("ClearFiber",0.,clrfiberHalfL,clrfiberZ,0.0*rad,twopi*rad);

  G4LogicalVolume*   logicClrfiber =
                                          new G4LogicalVolume(solidClrfiber,
                                                       FindMaterial("G4_AIR"),
                                                       "ClearFiber");

  new G4PVPlacement(new G4RotationMatrix(CLHEP::HepRotationY(-mppcTheta)),
                    G4ThreeVector(mppcOriginX,0.0,mppcOriginZ),
                    logicClrfiber,
                    "ClearFiber",
                    logicCouple,
                    false,
                    0); 

  //--------------------------------------------------
  // PhotonDet (Sensitive Detector)
  //--------------------------------------------------  

  // Physical Construction
  G4VSolid* solidPhotonDet;

  if ( mppcShape == "Square" )
    solidPhotonDet = new G4Box("PhotonDet",mppcHalfL,mppcHalfL,mppcZ);
  else
    solidPhotonDet =
                 new G4Tubs("PhotonDet",0.,mppcHalfL,mppcZ,0.0*rad,twopi*rad);

  G4LogicalVolume*   logicPhotonDet =
                                    new G4LogicalVolume(solidPhotonDet,
                                                        FindMaterial("G4_Al"),
                                                        "PhotonDet");

  new G4PVPlacement(0,
                    G4ThreeVector(0.0,0.0,0.0),
                    logicPhotonDet,
                    "PhotonDet",
                    logicClrfiber,
                    false,
                    0);

  // PhotonDet Surface Properties
  G4OpticalSurface* PhotonDetSurface = new G4OpticalSurface("PhotonDetSurface",
                                                       glisur,
                                                       ground,
                                                       dielectric_metal,
                                                       mppcPolish);

  G4MaterialPropertiesTable* PhotonDetSurfaceProperty =
                                               new G4MaterialPropertiesTable();

  G4double p_mppc[2] = {2.00*eV, 3.47*eV};
  G4double refl_mppc[2] = {mppcReflectivity,mppcReflectivity};
  G4double effi_mppc[2] = {1, 1};
 
  PhotonDetSurfaceProperty -> AddProperty("REFLECTIVITY",p_mppc,refl_mppc,2);
  PhotonDetSurfaceProperty -> AddProperty("EFFICIENCY",p_mppc,effi_mppc,2);

  PhotonDetSurface -> SetMaterialPropertiesTable(PhotonDetSurfaceProperty);
 
  new G4LogicalSkinSurface("PhotonDetSurface",logicPhotonDet,PhotonDetSurface); 

  if (!mppcSD) {
     G4String mppcSDName = "WLS/PhotonDet";
     mppcSD = new WLSPhotonDetSD(mppcSDName);

     G4SDManager* SDman = G4SDManager::GetSDMpointer();
     SDman->AddNewDetector(mppcSD);
  }

  // Setting the detector to be sensitive
  logicPhotonDet->SetSensitiveDetector(mppcSD);

}

void WLSDetectorConstruction::UpdateGeometry()
{
  if (!physiWorld) return;

  // clean-up previous geometry
  G4GeometryManager::GetInstance()->OpenGeometry();

  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  G4LogicalSkinSurface::CleanSurfaceTable();
  G4LogicalBorderSurface::CleanSurfaceTable();

  //define new one
  UpdateGeometryParameters();
 
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetector());

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();

  G4RegionStore::GetInstance()->UpdateMaterialList(physiWorld);
}

void WLSDetectorConstruction::UpdateGeometryParameters()
{
  wlsfiberRX  = XYRatio * wlsfiberRY;

  clad1RX = wlsfiberRX + 0.03*wlsfiberRX;
  clad1RY = wlsfiberRY + 0.03*wlsfiberRY;
  clad1Z  = wlsfiberZ;

  clad2RX = clad1RX + 0.03*wlsfiberRX;
  clad2RY = clad1RY + 0.03*wlsfiberRY;
  clad2Z  = wlsfiberZ;

  worldSizeX = clad2RX   + mppcDist + mppcHalfL + 1.*cm;
  worldSizeY = clad2RY   + mppcDist + mppcHalfL + 1.*cm;
  worldSizeZ = wlsfiberZ + mppcDist + mppcHalfL + 1.*cm;
 
  coupleRX = worldSizeX;
  coupleRY = worldSizeY;
  coupleZ  = (worldSizeZ - wlsfiberZ) / 2;
 
  clrfiberHalfL = mppcHalfL;
 
  mirrorRmax = clad2RY;
 
  coupleOrigin = wlsfiberOrigin + wlsfiberZ + coupleZ; 
  mirrorOrigin = wlsfiberOrigin - wlsfiberZ - mirrorZ; 
  mppcOriginX  = std::sin(mppcTheta) * (mppcDist + clrfiberZ);
  mppcOriginZ  = -coupleZ + std::cos(mppcTheta) * (mppcDist + clrfiberZ);
}

G4RotationMatrix 
            WLSDetectorConstruction::stringToRotationMatrix(G4String rotation)
{
  // We apply successive rotations OF THE OBJECT around the FIXED
  // axes of the parent's local coordinates; rotations are applied
  // left-to-right (rotation="r1,r2,r3" => r1 then r2 then r3).

  G4RotationMatrix rot;

  unsigned int place = 0;

  while (place < rotation.size()) {

        G4double angle;
        char* p;

        angle = strtod(rotation.substr(place+1).c_str(),&p) * deg;

        if (!p || (*p != (char)',' && *p != (char)'\0')) { 
           G4cerr << "Invalid rotation specification: " << 
                                                  rotation.c_str() << G4endl;
           return rot;
        }

        G4RotationMatrix thisRotation;

        switch(rotation.substr(place,1).c_str()[0]) {
              case 'X': case 'x':
                thisRotation = G4RotationMatrix(CLHEP::HepRotationX(angle));
                break;
              case 'Y': case 'y':
                thisRotation = G4RotationMatrix(CLHEP::HepRotationY(angle));
                break;
              case 'Z': case 'z':
                thisRotation = G4RotationMatrix(CLHEP::HepRotationZ(angle));
                break;
              default:
                G4cerr << " Invalid rotation specification: "
                       << rotation << G4endl;
                return rot;
        }

       rot = thisRotation * rot;
       place = rotation.find(',',place);
       if (place > rotation.size()) break;
       ++place;
  }

  return rot;

}

// Set the Geometry of the PhotonDet detector
// Pre:  shape must be either "Circle" and "Square"
void WLSDetectorConstruction::SetPhotonDetGeometry (G4String shape)
{
  if (shape == "Circle" || shape == "Square" ) mppcShape = shape;
}

// Set the number of claddings
// Pre: 0 <= num <= 2
void WLSDetectorConstruction::SetNumberOfCladding(G4int num)
{
  numOfCladLayers = num;
}

// Set the TOTAL length of the WLS fiber
void WLSDetectorConstruction::SetWLSLength (G4double length)
{
  wlsfiberZ = length;
}

// Set the Y radius of WLS fiber
void WLSDetectorConstruction::SetWLSRadius (G4double radius)
{
  wlsfiberRY = radius;
}

// Set the Y radius of Cladding 1
void WLSDetectorConstruction::SetClad1Radius (G4double radius)
{
  clad1RY = radius;
}

// Set the Y radius of Cladding 2
void WLSDetectorConstruction::SetClad2Radius (G4double radius)
{
  clad2RY = radius;
}

// Set the half length of the PhotonDet detector
// The half length will be the radius if PhotonDet is circular
void WLSDetectorConstruction::SetPhotonDetHalfLength(G4double halfL)
{
  mppcHalfL = halfL;
}

// Set the distance between fiber end and PhotonDet
void WLSDetectorConstruction::SetGap (G4double gap) { mppcDist = gap; }

// Set the Aligment of PhotonDet with respect to the z axis
// If theta is 0 deg, then the detector is perfectly aligned
// PhotonDet will be deviated by theta from z axis
// facing towards the center of the fiber
void WLSDetectorConstruction::SetPhotonDetAlignment(G4double theta)
{
  mppcTheta = theta;
}

// Set the Surface Roughness between Cladding 1 and WLS fiber
// Pre: 0 < roughness <= 1
void WLSDetectorConstruction::SetSurfaceRoughness(G4double roughness)
{
  surfaceRoughness = roughness;
}

// Set the Polish of the mirror, polish of 1 is a perfect mirror surface
// Pre: 0 < polish <= 1
void WLSDetectorConstruction::SetMirrorPolish(G4double polish)
{
  mirrorPolish = polish;
}

// Set the Reflectivity of the mirror, reflectivity of 1 is a perfect mirror
// Pre: 0 < reflectivity <= 1
void WLSDetectorConstruction::SetMirrorReflectivity(G4double reflectivity)
{
  mirrorReflectivity = reflectivity;
}

// Set the Polish of the PhotonDet, polish of 1 is a perfect mirror surface
// Pre: 0 < polish <= 1
void WLSDetectorConstruction::SetPhotonDetPolish(G4double polish)
{
  mppcPolish = polish;
}

// Set the Reflectivity of the PhotonDet, reflectivity of 1 is a perfect mirror
// Pre: 0 < reflectivity <= 1
void WLSDetectorConstruction::SetPhotonDetReflectivity(G4double reflectivity)
{
  mppcReflectivity = reflectivity;
}

// Toggle to place the mirror or not at one end (-z end) of the fiber
// True means place the mirror, false means otherwise
void WLSDetectorConstruction::SetMirror(G4bool flag) { mirrorToggle = flag; }

// Set the ratio of the x and y radius of the ellipse (x/y)
// a ratio of 1 would produce a circle
void WLSDetectorConstruction::SetXYRatio(G4double r) { XYRatio = r; }

// Set the length of the scintillator bar
void WLSDetectorConstruction::SetBarLength (G4double length)
{
  barLength = length;
}

// Set the side of the scintillator bar
void WLSDetectorConstruction::SetBarBase (G4double side)
{
  barBase = side;
}

// Set the radius of the fiber hole
void WLSDetectorConstruction::SetHoleRadius (G4double radius)
{
  holeRadius = radius;
}

// Set thickness of the coating on the bars
void WLSDetectorConstruction::SetCoatingThickness (G4double thick)
{
  coatingThickness = thick;
}

// Set inner radius of the corner bar coating
void WLSDetectorConstruction::SetCoatingRadius (G4double radius)
{
  coatingRadius = radius;
}

G4double WLSDetectorConstruction::GetWLSFiberLength() { return wlsfiberZ; }

G4double WLSDetectorConstruction::GetBarLength() { return barLength; }
G4double WLSDetectorConstruction::GetBarBase() { return barBase; }
G4double WLSDetectorConstruction::GetHoleRadius() { return holeRadius; }
G4double WLSDetectorConstruction::GetHoleLength() { return holeLength; }
G4double WLSDetectorConstruction::GetFiberRadius() { return GetWLSFiberRMax(); }
G4double WLSDetectorConstruction::GetCoatingThickness() 
                                                  { return coatingThickness; }
G4double WLSDetectorConstruction::GetCoatingRadius() { return coatingRadius; }

G4double WLSDetectorConstruction::GetWLSFiberEnd()
{
  return wlsfiberOrigin + wlsfiberZ;
}

G4double WLSDetectorConstruction::GetWLSFiberRMax()
{
  if (numOfCladLayers == 2) return clad2RY;
  if (numOfCladLayers == 1) return clad1RY;
  return wlsfiberRY;
}

G4double WLSDetectorConstruction::GetSurfaceRoughness()
{
  return surfaceRoughness;
}

// Return True if the fiber construction is ideal
G4bool WLSDetectorConstruction::IsPerfectFiber()
{ 
  return     surfaceRoughness == 1. && XYRatio == 1.
             && (!mirrorToggle    || 
             (mirrorPolish    == 1. && mirrorReflectivity == 1.));
}

G4Material* WLSDetectorConstruction::FindMaterial(G4String name) {
    G4Material* material = G4Material::GetMaterial(name,true);
    return material;
}

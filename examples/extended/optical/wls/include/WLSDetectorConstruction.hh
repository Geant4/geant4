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
/// \file optical/wls/include/WLSDetectorConstruction.hh
/// \brief Definition of the WLSDetectorConstruction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef WLSDetectorConstruction_h
#define WLSDetectorConstruction_h 1

#include "globals.hh"
#include "G4Cache.hh"
#include "G4RotationMatrix.hh"
#include "G4VUserDetectorConstruction.hh"
#include <CLHEP/Units/SystemOfUnits.h>

class WLSMaterials;
class WLSDetectorMessenger;
class WLSPhotonDetSD;

class G4Box;
class G4EllipticalTube;
class G4LogicalVolume;
class G4Material;
class G4Tubs;
class G4VisAttributes;
class G4VPhysicalVolume;

class WLSDetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  WLSDetectorConstruction();
  ~WLSDetectorConstruction() override;

  G4VPhysicalVolume* Construct() override;
  G4VPhysicalVolume* ConstructDetector();

  void ConstructSDandField() override;

  // Set Material Commands for World and WLSfiber
  void SetWorldMaterial(G4String);
  void SetWLSFiberMaterial(G4String);
  void SetCoupleMaterial(G4String);

  void SetPhotonDetGeometry(G4String);
  void SetNumberOfCladding(G4int);  // Maximum 2 claddings

  void SetWLSLength(G4double);  // Total length of WLS fiber
  void SetWLSRadius(G4double);
  void SetClad1Radius(G4double);
  void SetClad2Radius(G4double);
  void SetPhotonDetHalfLength(G4double);
  void SetGap(G4double);
  void SetPhotonDetAlignment(G4double);
  // Set the ratio of x and y (x/y) radius of the ellipse
  void SetXYRatio(G4double);
  // Set the Roughness in between each layer
  void SetSurfaceRoughness(G4double);
  // Set the reflectivity of the mirror
  void SetMirrorReflectivity(G4double);
  // Set the polish of the mirror
  void SetMirrorPolish(G4double);
  // Set the reflectivity of the mirror
  void SetPhotonDetReflectivity(G4double);
  // Set the polish of the mirror
  void SetPhotonDetPolish(G4double);

  void SetMirror(G4bool);

  void SetBarLength(G4double);
  void SetBarBase(G4double);
  void SetHoleRadius(G4double);
  void SetCoatingThickness(G4double);
  void SetCoatingRadius(G4double);

  G4double GetWLSFiberLength();
  G4double GetWLSFiberEnd();
  G4double GetWLSFiberRMax();
  G4double GetSurfaceRoughness();
  G4bool IsPerfectFiber();

  G4double GetBarLength();
  G4double GetBarBase();
  G4double GetHoleRadius();
  G4double GetHoleLength();
  G4double GetFiberRadius();

  G4double GetCoatingThickness();
  G4double GetCoatingRadius();

  G4Material* FindMaterial(G4String);

 private:
  std::vector<G4VisAttributes*> fVisAttributes;

  WLSMaterials* fMaterials = nullptr;

  G4LogicalVolume* fLogicHole = nullptr;
  G4LogicalVolume* fLogicWorld = nullptr;

  G4VPhysicalVolume* fPhysiWorld = nullptr;
  G4VPhysicalVolume* fPhysiHole = nullptr;

  G4double fWorldSizeX = -1.;
  G4double fWorldSizeY = -1.;
  G4double fWorldSizeZ = -1.;

  G4double fWLSfiberRX = -1.;
  G4double fWLSfiberRY = 0.5 * CLHEP::mm;
  G4double fWLSfiberZ = 1. * CLHEP::m;

  G4double fClad1RX = -1.;
  G4double fClad1RY = -1.;
  G4double fClad1Z = -1.;

  G4double fClad2RX = -1.;
  G4double fClad2RY = -1.;
  G4double fClad2Z = -1.;

  G4double fClrfiberHalfL = -1.;
  G4double fClrfiberZ = -1.;

  G4double fCoupleRX = -1.;
  G4double fCoupleRY = -1.;
  G4double fCoupleZ = -1.;

  G4double fMirrorRmax = -1.;
  G4double fMirrorZ = 0.1 * CLHEP::mm;
  G4bool fMirrorToggle = true;

  G4String fMPPCShape = "Circle";
  G4double fMPPCHalfL = -1.;
  G4double fMPPCZ = 0.05 * CLHEP::mm;
  G4double fMPPCDist = 0.;
  G4double fMPPCTheta = 0;

  G4double fWLSfiberOrigin = 0.;
  G4double fCoupleOrigin = 0.;
  G4double fMirrorOrigin = 0.;
  G4double fMPPCOriginX = 0.;
  G4double fMPPCOriginZ = 0.;

  G4int fNumOfCladLayers = 0;

  G4double fMirrorPolish = 1.;
  G4double fMirrorReflectivity = 1.;
  G4double fMPPCPolish = 1.;
  G4double fMPPCReflectivity = 0.;
  G4double fExtrusionPolish = 1.;
  G4double fExtrusionReflectivity = 1.;
  G4double fSurfaceRoughness = 1.;
  G4double fXYRatio = 1.;

  G4double fBarLength = 1. * CLHEP::m;
  G4double fBarBase = 9.6 * CLHEP::mm;
  G4double fHoleRadius = 0.9 * CLHEP::mm;
  G4double fHoleLength = -1.;
  G4double fCoatingThickness = 0.25 * CLHEP::mm;
  G4double fCoatingRadius = 1.875 * CLHEP::mm;

  void UpdateGeometryParameters();

  WLSDetectorMessenger* fDetectorMessenger = nullptr;
  G4Cache<WLSPhotonDetSD*> fmppcSD;
};

#endif

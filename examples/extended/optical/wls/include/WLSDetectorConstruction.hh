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
  ~WLSDetectorConstruction();

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

  WLSMaterials* fMaterials;

  G4LogicalVolume* fLogicHole;
  G4LogicalVolume* fLogicWorld;

  G4VPhysicalVolume* fPhysiWorld;
  G4VPhysicalVolume* fPhysiHole;

  G4double fWorldSizeX;
  G4double fWorldSizeY;
  G4double fWorldSizeZ;

  G4double fWLSfiberRX;
  G4double fWLSfiberRY;
  G4double fWLSfiberZ;

  G4double fClad1RX;
  G4double fClad1RY;
  G4double fClad1Z;

  G4double fClad2RX;
  G4double fClad2RY;
  G4double fClad2Z;

  G4double fClrfiberHalfL;
  G4double fClrfiberZ;

  G4double fCoupleRX;
  G4double fCoupleRY;
  G4double fCoupleZ;

  G4double fMirrorRmax;
  G4double fMirrorZ;
  G4bool fMirrorToggle;

  G4String fMPPCShape;
  G4double fMPPCHalfL;
  G4double fMPPCZ;
  G4double fMPPCDist;
  G4double fMPPCTheta;

  G4double fWLSfiberOrigin;
  G4double fCoupleOrigin;
  G4double fMirrorOrigin;
  G4double fMPPCOriginX;
  G4double fMPPCOriginZ;

  G4int fNumOfCladLayers;

  G4double fMirrorPolish;
  G4double fMirrorReflectivity;
  G4double fMPPCPolish;
  G4double fMPPCReflectivity;
  G4double fExtrusionPolish;
  G4double fExtrusionReflectivity;
  G4double fSurfaceRoughness;
  G4double fXYRatio;

  G4double fBarLength;
  G4double fBarBase;
  G4double fHoleRadius;
  G4double fHoleLength;
  G4double fCoatingThickness;
  G4double fCoatingRadius;

  void UpdateGeometryParameters();

  WLSDetectorMessenger* fDetectorMessenger;
  G4Cache<WLSPhotonDetSD*> fmppcSD;
};

#endif

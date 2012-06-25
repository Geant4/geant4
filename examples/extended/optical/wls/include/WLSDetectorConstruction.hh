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
/// \file optical/wls/include/WLSDetectorConstruction.hh
/// \brief Definition of the WLSDetectorConstruction class
//
//
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef WLSDetectorConstruction_h
#define WLSDetectorConstruction_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4RotationMatrix.hh"

class G4Box;
class G4Tubs;
class G4EllipticalTube;

class G4LogicalVolume;
class G4VPhysicalVolume;

class WLSMaterials;
class G4Material;

class WLSDetectorMessenger;

class WLSPhotonDetSD;

#include "G4VUserDetectorConstruction.hh"

class WLSDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    WLSDetectorConstruction();
    ~WLSDetectorConstruction();

    G4VPhysicalVolume* Construct();
    G4VPhysicalVolume* ConstructDetector();

    void UpdateGeometry();
 
    // Set Material Commands for World and WLSfiber
    void SetWorldMaterial         (G4String);
    void SetWLSFiberMaterial      (G4String);
    void SetCoupleMaterial        (G4String);

    void SetPhotonDetGeometry     (G4String);
    void SetNumberOfCladding      (G4int);        // Maximum 2 claddings

    void SetWLSLength             (G4double);     // Total length of WLS fiber
    void SetWLSRadius             (G4double);
    void SetClad1Radius           (G4double);
    void SetClad2Radius           (G4double);
    void SetPhotonDetHalfLength   (G4double);
    void SetGap                   (G4double);
    void SetPhotonDetAlignment    (G4double);
    // Set the ratio of x and y (x/y) radius of the ellipse
    void SetXYRatio               (G4double);
    // Set the Roughness in between each layer
    void SetSurfaceRoughness      (G4double);
    // Set the reflectivity of the mirror
    void SetMirrorReflectivity    (G4double);
    // Set the polish of the mirror
    void SetMirrorPolish          (G4double);
    // Set the reflectivity of the mirror
    void SetPhotonDetReflectivity (G4double);
    // Set the polish of the mirror
    void SetPhotonDetPolish       (G4double);

    void SetMirror                (G4bool);

    void SetBarLength             (G4double);
    void SetBarBase               (G4double);
    void SetHoleRadius            (G4double);
    void SetCoatingThickness      (G4double);
    void SetCoatingRadius         (G4double);

    G4double GetWLSFiberLength();
    G4double GetWLSFiberEnd();
    G4double GetWLSFiberRMax();
    G4double GetSurfaceRoughness();
    G4bool   IsPerfectFiber();

    G4double GetBarLength();
    G4double GetBarBase();
    G4double GetHoleRadius();
    G4double GetHoleLength();
    G4double GetFiberRadius();

    G4double GetCoatingThickness();
    G4double GetCoatingRadius();
 
    // stringToRotationMatrix() converts a string "X90,Y45" into a
    // G4RotationMatrix.
    // This is an active rotation, in that the object is first rotated
    // around the parent's X axis by 90 degrees, then the object is
    // further rotated around the parent's Y axis by 45 degrees.
    // The return value points to a G4RotationMatrix on the heap, so
    // it is persistent. Angles are in degrees, can have decimals,
    // and can be negative. Axes are X, Y, Z.

    static G4RotationMatrix stringToRotationMatrix(G4String rotation);

    G4Material* FindMaterial(G4String);

  private:

    WLSMaterials* materials;

    G4LogicalVolume* logicHole;
    G4LogicalVolume* logicWorld;

    G4VPhysicalVolume* physiWorld;
    G4VPhysicalVolume* physiHole;
 
    G4double           worldSizeX;
    G4double           worldSizeY;
    G4double           worldSizeZ;

    G4double           wlsfiberRX;
    G4double           wlsfiberRY;
    G4double           wlsfiberZ;

    G4double           clad1RX;
    G4double           clad1RY;
    G4double           clad1Z;

    G4double           clad2RX;
    G4double           clad2RY;
    G4double           clad2Z;

    G4double           clrfiberHalfL;
    G4double           clrfiberZ;

    G4double           coupleRX;
    G4double           coupleRY;
    G4double           coupleZ;
 
    G4double           mirrorRmax;
    G4double           mirrorZ;
    G4bool             mirrorToggle;
 
    G4String           mppcShape;
    G4double           mppcHalfL;
    G4double           mppcZ;
    G4double           mppcDist;
    G4double           mppcTheta;

    G4double wlsfiberOrigin;
    G4double coupleOrigin;
    G4double mirrorOrigin;
    G4double mppcOriginX;
    G4double mppcOriginZ;
 
    G4int numOfCladLayers;

    G4double mirrorPolish;
    G4double mirrorReflectivity;
    G4double mppcPolish;
    G4double mppcReflectivity;
    G4double extrusionPolish;
    G4double extrusionReflectivity;
    G4double surfaceRoughness;
    G4double XYRatio;

    G4double barLength;
    G4double barBase;
    G4double holeRadius;
    G4double holeLength;
    G4double coatingThickness;
    G4double coatingRadius;

    static WLSPhotonDetSD* mppcSD;

  private:

     void ConstructFiber();

     void UpdateGeometryParameters();

     WLSDetectorMessenger* detectorMessenger;

};

#endif

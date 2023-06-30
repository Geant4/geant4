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

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4Cache.hh"
#include "G4LogicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4Ellipsoid;
class G4Sphere;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;
class G4GlobalMagFieldMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  ~DetectorConstruction() override;

  void SetAbsorberMaterial(const G4String&);
  void SetAbsorberThickness(G4double);
  void SetAbsorberSizeYZ(G4double);

  void SetAbsorberXpos(G4double);

  void SetWorldMaterial(const G4String&);
  void SetWorldSizeX(G4double);
  void SetWorldSizeYZ(G4double);
  void SetPhantomType(G4int value);

  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;

  void PrintGeomParameters();

private:
  void DefineMaterials();
  void ComputeGeomParameters();
  void ChangeGeometry();
  void Construct_Phantom1();
  void Construct_Phantom2();
  void Construct_Phantom3();

  G4int phantom_type;

  G4Material* fAbsorberMaterial;
  G4double fAbsorberThickness;
  G4double fAbsorberSizeYZ;

  G4double fXposAbs;
  G4double fXstartAbs, fXendAbs;

  G4Material* fWorldMaterial;
  G4double fWorldSizeX;
  G4double fWorldSizeYZ;

  G4Box* fSolidWorld;
  G4LogicalVolume* fLogicWorld;
  G4VPhysicalVolume* fPhysiWorld;

  G4Box* fSolidAbsorber;
  G4LogicalVolume* fLogicAbsorber;
  G4VPhysicalVolume* fPhysiAbsorber;

  G4Material* material1;
  G4Material* material2;
  G4Material* material3;
  G4Material* material4;
  G4Material* material5;
  G4Material* material6;
  G4Material* material_GDP;

  G4Ellipsoid* ellipse1;
  G4LogicalVolume* logicEllipse1;
  G4VPhysicalVolume* physiEllipse1;

  G4Ellipsoid* ellipse2;
  G4LogicalVolume* logicEllipse2;
  G4VPhysicalVolume* physiEllipse2;

  G4Ellipsoid* ellipse3;
  G4LogicalVolume* logicEllipse3;
  G4VPhysicalVolume* physiEllipse3;

  G4Ellipsoid* ellipse4;
  G4LogicalVolume* logicEllipse4;
  G4VPhysicalVolume* physiEllipse4;

  G4Ellipsoid* ellipse5;
  G4LogicalVolume* logicEllipse5;
  G4VPhysicalVolume* physiEllipse5;

  G4Ellipsoid* ellipse6;
  G4LogicalVolume* logicEllipse6;
  G4VPhysicalVolume* physiEllipse6;

  // Ge-doped glow discharge polymer (GDP)
  G4Sphere* solid_GDP;
  G4LogicalVolume* logic_GDP;
  G4VPhysicalVolume* physi_GDP;

  DetectorMessenger* fDetectorMessenger;
  G4Cache<G4GlobalMagFieldMessenger*> fFieldMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

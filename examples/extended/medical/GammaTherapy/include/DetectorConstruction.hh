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
#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

// -------------------------------------------------------------
//      GEANT4  test  IBREM
//
// Authors: V.Grichine, V.Ivanchenko
//
// Modified:
//
// 18-02-03 V.Ivanchenko create
//
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class CheckVolumeSD;
class PhantomSD;
class TargetSD;
class DetectorMessenger;
class G4LogicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  DetectorConstruction();
  ~DetectorConstruction();

public:

  G4VPhysicalVolume* Construct();
  G4double GetGeneratorPosZ() const {return fGeneratorPosZ;};

  void setGap(G4double val) {fDelta = val;};
  void setTarget1Z(G4double val) {fTarget1Z = val;};
  void setTarget2Z(G4double val) {fTarget2Z = val;};
  void setMylarZ(G4double val) {fMylarVolumeZ = val;};
  void setCheckShiftZ(G4double val) {fCheckShiftZ = val;};
  void setAbsorberZ(G4double val) {fAbsorberZ = val;};
  void setAbsorberShiftZ(G4double val) {fAbsorberShiftZ = val;};

  void setTarget1Material(const G4String& m);
  void setTarget2Material(const G4String& m);
  void UpdateGeometry();

private:  // methods

  void Initialise();
  void DefineMaterials();
  void InitialiseGeometryParameters();
  G4VPhysicalVolume* ConstructVolumes();

private:  // fields

  CheckVolumeSD* checkSD;
  PhantomSD* calorimeterSD;
  TargetSD* targetSD;
  G4double fWorldXY, fWorldZ;
  G4double fDelta;
  G4double fGeneratorPosZ;

  G4double fTargetRadius, fTarget1Z, fTarget1PosZ;
  G4double fTarget2Z, fTarget2PosZ;

  G4double fGasVolumeRadius,fGasVolumeZ, fGasVolumePosZ, fAirZ, fMylarVolumeZ, fMylarPosZ;
  G4double fCheckVolumeRadius, fCheckVolumeZ, fCheckShiftZ, fCheckVolumePosZ;
  G4double fTargetVolumeZ, fTargetVolumePosZ;

  G4double fPhantomRadius, fPhantomZ, fPhantomPosZ;
  G4double fAbsorberRadius, fAbsorberZ, fAbsorberShiftZ, fAbsorberPosZ;
  G4double fDistanceVacuumTarget, fWindowZ, fWindowPosZ;

  G4Material* fWorldMaterial;

  G4Material* fTarget1Material;
  G4Material* fTarget2Material;
  G4Material* fMylar;
  G4Material* fWindowMaterial;

  G4Material* fLightMaterial;
  G4Material* fAbsorberMaterial;

  G4LogicalVolume* logicTarget1;
  G4LogicalVolume* logicTarget2;

  DetectorMessenger*  dMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

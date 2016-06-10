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
// $Id: DetectorConstruction.hh 67994 2013-03-13 11:05:39Z gcosmo $
//
/// \file medical/GammaTherapy/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
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
  virtual ~DetectorConstruction();

  G4VPhysicalVolume* Construct();

  void SetTarget1Material(const G4String& m);
  void SetTarget2Material(const G4String& m);
  void UpdateGeometry();

  inline G4double GetGeneratorPosZ() const    { return fGeneratorPosZ; };

  inline void SetGap(G4double val)            { fDelta = val; 
                                                UpdateGeometry(); };
  inline void SetTarget1Z(G4double val)       { fTarget1Z = val; 
                                                UpdateGeometry(); };
  inline void SetTarget2Z(G4double val)       { fTarget2Z = val; 
                                                UpdateGeometry(); };
  inline void SetMylarZ(G4double val)         { fMylarVolumeZ = val; 
                                                UpdateGeometry(); };
  inline void SetCheckShiftZ(G4double val)    { fCheckShiftZ = val; 
                                                UpdateGeometry(); };
  inline void SetAbsorberZ(G4double val)      { fAbsorberZ = val; 
                                                UpdateGeometry(); };
  inline void SetAbsorberShiftZ(G4double val) { fAbsorberShiftZ = val; 
                                                UpdateGeometry(); };

private:  

  void InitialiseGeometryParameters();

  DetectorConstruction & operator=(const DetectorConstruction &right);
  DetectorConstruction(const DetectorConstruction&);

  CheckVolumeSD* fCheckSD;
  PhantomSD*     fPhantomSD;
  TargetSD*      fTargetSD;

  G4double fWorldXY, fWorldZ;
  G4double fDelta;
  G4double fGeneratorPosZ;

  G4double fTargetRadius, fTarget1Z, fTarget1PosZ;
  G4double fTarget2Z, fTarget2PosZ;

  G4double fGasVolumeRadius, fGasVolumeZ, fGasVolumePosZ;
  G4double fAirZ, fMylarVolumeZ, fMylarPosZ;
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

  G4LogicalVolume* fLogicTarget1;
  G4LogicalVolume* fLogicTarget2;

  DetectorMessenger* fMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

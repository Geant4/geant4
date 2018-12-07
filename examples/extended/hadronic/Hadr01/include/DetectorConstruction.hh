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
/// \file hadronic/Hadr01/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
//
/////////////////////////////////////////////////////////////////////////
//
// DetectorConstruction
//
// Created: 31.01.2003 V.Ivanchenko
//
// Modified:
// 04.06.2006 Adoptation of Hadr01 (V.Ivanchenko)
//
////////////////////////////////////////////////////////////////////////
// 

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Material.hh"

class CheckVolumeSD;
class TargetSD;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  DetectorConstruction();
  virtual ~DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();

  void SetWorldMaterial(const G4String&);
  void SetTargetMaterial(const G4String&);

  void SetTargetRadius(G4double val);

private:

  void ComputeGeomParameters();

  DetectorConstruction & operator=(const DetectorConstruction &right);
  DetectorConstruction(const DetectorConstruction&);

  G4Material*  fTargetMaterial;
  G4Material*  fWorldMaterial;

  G4Tubs* fSolidW;
  G4Tubs* fSolidA;
  G4Tubs* fSolidC;

  G4LogicalVolume* fLogicTarget;
  G4LogicalVolume* fLogicCheck;
  G4LogicalVolume* fLogicWorld;

  G4VPhysicalVolume* fPhysWorld;

  DetectorMessenger* fDetectorMessenger;

  G4double fRadius;
  G4double fCheckR;
  G4double fWorldR;
  G4double fTargetZ; 
  G4double fCheckZ;
  G4double fWorldZ;
  G4double fSliceZ;

  G4int fSlices;
  G4bool fInitialized;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#endif


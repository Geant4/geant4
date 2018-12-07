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
/// \file electromagnetic/TestEm9/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
//
//
/////////////////////////////////////////////////////////////////////////
//
// TestEm9: Cryslal calorimetry
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"

class G4Tubs;
class G4LogicalVolume;
class G4UniformMagField;
class DetectorMessenger;
class G4Region;
class G4ProductionCuts;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  DetectorConstruction();
  virtual ~DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();

  void SetEcalMaterial(const G4String&);
  void SetAbsMaterial(const G4String&);
  void SetEcalLength(G4double val);
  void SetEcalWidth(G4double val);
  void SetVertexLength(G4double val);
  void SetPadLength(G4double val);
  void SetPadWidth(G4double val);
  void SetAbsLength(G4double val);

  void UpdateGeometry();

  inline G4double GetWorldSizeZ()  {return fWorldZ;}

private:

  void DefineMaterials();
  G4VPhysicalVolume* ConstructVolumes();

  DetectorConstruction & operator=(const DetectorConstruction &right);
  DetectorConstruction(const DetectorConstruction&);

  G4double fEcalLength;
  G4double fEcalWidth;
  G4double fVertexLength;
  G4double fPadLength;
  G4double fPadWidth;
  G4double fAbsLength;
  G4double fWorldZ;

  G4Material* fCalMaterial;
  G4Material* fVertMaterial;
  G4Material* fAbsMaterial;
  G4Material* fWorldMaterial;
  G4Material* fYorkMaterial;

  G4LogicalVolume* fLogicWorld;
  G4LogicalVolume* fLogicCal;
  G4LogicalVolume* fLogicA1;
  G4LogicalVolume* fLogicA2;
  G4LogicalVolume* fLogicA3;
  G4LogicalVolume* fLogicA4;

  G4Region* fVertexRegion;
  G4Region* fMuonRegion;
  
  G4ProductionCuts* fVertexDetectorCuts;
  G4ProductionCuts* fMuonDetectorCuts;

  DetectorMessenger* fDetectorMessenger;  //pointer to the Messenger
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#endif


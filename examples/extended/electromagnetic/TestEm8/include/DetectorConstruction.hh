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
/// \file electromagnetic/TestEm8/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
//
/////////////////////////////////////////////////////////////////////////
//
// TestEm8: Gaseous detector
//
// Created: 31.08.2010 V.Ivanchenko on base of V.Grichine code
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Material;
class G4Region;
class DetectorMessenger;
class TargetSD;
class G4Tubs;
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4ProductionCuts;

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  explicit DetectorConstruction();
  ~DetectorConstruction() override;

  G4VPhysicalVolume *Construct() override;
  void ConstructSDandField() override;

  void SetGasMaterial(const G4String &);
  void SetContainerMaterial(const G4String &);
  void SetWorldMaterial(const G4String &);

  void SetGasThickness(G4double);
  void SetGasRadius(G4double);
  void SetContainerThickness(G4double);

  void SetPairEnergy(G4double);

  DetectorConstruction &operator = (const DetectorConstruction &right) = delete;
  DetectorConstruction(const DetectorConstruction &) = delete;

  inline G4VPhysicalVolume *GetWorldPhysVol() const { return fPhysWorld; }
  inline void SetMaxChargedStep(G4double x) { fMaxStep = x; }
  inline G4double GetMaxChargedStep() const { return fMaxStep; }

private:
  void DefineMaterials();
  void ChangeGeometry();

  G4Material *fGasMat;
  G4Material *fWindowMat;
  G4Material *fWorldMaterial;

  G4double fGasThickness;
  G4double fGasRadius;
  G4double fMaxStep;
  G4double fWindowThick;


  G4Tubs *fSolidWorld = nullptr;
  G4Tubs *fSolidContainer = nullptr;
  G4Tubs *fSolidDetector = nullptr;

  G4VPhysicalVolume *fPhysWorld = nullptr;
  G4LogicalVolume *fLogicWorld = nullptr;
  G4LogicalVolume *fLogicContainer = nullptr;
  G4LogicalVolume *fLogicDetector = nullptr;

  DetectorMessenger *fDetectorMessenger;
  G4ProductionCuts *fGasDetectorCuts;
  G4Region *fRegGasDet = nullptr;
};

#endif

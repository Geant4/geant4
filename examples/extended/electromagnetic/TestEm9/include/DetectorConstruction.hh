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
// $Id: DetectorConstruction.hh,v 1.7 2008-04-07 18:09:05 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

public:

  void SetEcalMaterial(const G4String&);
  void SetAbsMaterial(const G4String&);
  void SetEcalLength(G4double val);
  void SetEcalWidth(G4double val);
  void SetVertexLength(G4double val);
  void SetPadLength(G4double val);
  void SetPadWidth(G4double val);
  void SetAbsLength(G4double val);

  G4VPhysicalVolume* Construct();

  void UpdateGeometry();

  G4double GetWorldSizeZ()  {return worldZ;}

private:

  void DefineMaterials();
  G4VPhysicalVolume* ConstructVolumes();

  DetectorConstruction & operator=(const DetectorConstruction &right);
  DetectorConstruction(const DetectorConstruction&);

  G4double ecalLength;
  G4double ecalWidth;
  G4double vertexLength;
  G4double padLength;
  G4double padWidth;
  G4double absLength;
  G4double worldZ;

  G4Material* calMaterial;
  G4Material* vertMaterial;
  G4Material* absMaterial;
  G4Material* worldMaterial;
  G4Material* yorkMaterial;

  G4LogicalVolume* logicWorld;
  G4LogicalVolume* logicECal;
  G4LogicalVolume* logicCal;
  G4LogicalVolume* logicA1;
  G4LogicalVolume* logicA2;
  G4LogicalVolume* logicA3;
  G4LogicalVolume* logicA4;
  G4LogicalVolume* logicYV;
  G4LogicalVolume* logicY;
  G4LogicalVolume* logicVV;
  G4LogicalVolume* logicVD;
  G4LogicalVolume* logicV;
  
  G4Region*   vertexRegion;
  G4Region*   muonRegion;
  G4ProductionCuts* vertexDetectorCuts;
  G4ProductionCuts* muonDetectorCuts;

  DetectorMessenger* detectorMessenger;  //pointer to the Messenger
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#endif


//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: DetectorConstruction.hh,v 1.3 2004/05/27 13:43:17 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  DetectorConstruction();
  ~DetectorConstruction();

public:

  void SetEcalMaterial(const G4String&);
  void SetAbsMaterial(const G4String&);
  void SetEcalLength (G4double val)   {ecalLength = val;};
  void SetEcalWidth  (G4double val)   {ecalWidth = val;};
  void SetVertexLength (G4double val) {vertexLength = val;};
  void SetPadLength  (G4double val)   {padLength = val;};
  void SetPadWidth  (G4double val)    {padWidth = val;};
  void SetAbsLength(G4double val)     {absLength = val;};

  G4VPhysicalVolume* Construct();

  void UpdateGeometry();

  G4double GetWorldSizeZ()            {return worldZ;}

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

  G4LogicalVolume* logicC;
  G4LogicalVolume* logicA1;
  G4LogicalVolume* logicA2;
  G4LogicalVolume* logicA3;
  G4LogicalVolume* logicA4;
  
  G4Region*   vertexRegion;
  G4Region*   muonRegion;

  DetectorMessenger* detectorMessenger;  //pointer to the Messenger

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#endif


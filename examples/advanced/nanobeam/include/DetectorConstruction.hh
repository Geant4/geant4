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
// -------------------------------------------------------------------
// $Id: DetectorConstruction.hh,v 1.2 2008-01-25 20:49:24 sincerti Exp $
// -------------------------------------------------------------------

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EqMagElectricField.hh"
#include "G4UserLimits.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ClassicalRK4.hh"
#include "G4ChordFinder.hh"
#include "G4PropagatorInField.hh"

#include "DetectorMessenger.hh"
#include "TabulatedField3D.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  ~DetectorConstruction();

  G4VPhysicalVolume* Construct();
  void PrintDetectorParameters();
 
  void SetG1 (G4float);
  void SetG2 (G4float);
  void SetG3 (G4float);
  void SetG4 (G4float);
  void UpdateGeometry();

  void SetModel (G4int);
  
  G4int GetCoef();
  void SetCoef();
  
  G4int GetProfile(){return profile;};
  void SetProfile(G4int myProfile);
  
  G4int GetGrid(){return grid;};
  void SetGrid(G4int myGrid);

  G4float G1, G2, G3, G4;
  G4int model, coef, profile, grid;

private:
   
  void DefineMaterials();
  G4VPhysicalVolume* ConstructVolumes();     

  G4VPhysicalVolume* physiWorld;
  G4LogicalVolume*   logicWorld;  
  G4Box*             solidWorld;
  
  G4VPhysicalVolume* physiVol;
  G4LogicalVolume*   logicVol;  
  G4Box*             solidVol;
  
  G4VPhysicalVolume* physiGridVol;
  G4LogicalVolume*   logicGridVol;  
  G4Box*	     solidGridVol;

  G4VPhysicalVolume* physiGridVol_Hole;
  G4LogicalVolume*   logicGridVol_Hole;  
  G4Box*	     solidGridVol_Hole;

  G4VPhysicalVolume* physiControlVol_GridShadow;
  G4LogicalVolume*   logicControlVol_GridShadow;  
  G4Box*             solidControlVol_GridShadow;

  G4Material*        defaultMaterial;
  G4Material*        waterMaterial;
  G4Material*        gridMaterial;

  G4bool gradientsInitialized;
  DetectorMessenger* detectorMessenger;

};
#endif



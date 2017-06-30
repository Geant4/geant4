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
// Please cite the following paper if you use this software
// Nucl.Instrum.Meth.B260:20-27, 2007

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4TransportationManager.hh"
#include "G4UserLimits.hh"
#include "G4ClassicalRK4.hh"
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

  void ConstructSDandField();
     
  void PrintDetectorParameters();
 
  void SetG1 (G4float);
  void SetG2 (G4float);
  void SetG3 (G4float);
  void SetG4 (G4float);
  
  G4float GetG1 () {return fG1;};
  G4float GetG2 () {return fG2;};
  G4float GetG3 () {return fG3;};
  G4float GetG4 () {return fG4;};
  
  G4int GetModel(){return fModel;};
  void  SetModel (G4int);
  
  G4int GetCoef();
  void  SetCoef(G4int val);
  
  G4int GetProfile(){return fProfile;};
  void  SetProfile(G4int myProfile);
  
  G4int GetGrid(){return fGrid;};
  void  SetGrid(G4int myGrid);

  G4LogicalVolume* GetLogicalWorld() {return fLogicWorld;};
  G4LogicalVolume* GetLogicalVol() {return fLogicVol;};
  G4LogicalVolume* GetLogicalGrid() {return fLogicControlVol_GridShadow;};

  G4float fG1, fG2, fG3, fG4;
  G4int fModel, fCoef, fProfile, fGrid;

private:
   
  void DefineMaterials();
  G4VPhysicalVolume* ConstructVolumes();     

  G4VPhysicalVolume* fPhysiWorld;
  G4LogicalVolume*   fLogicWorld;  
  G4Box*             fSolidWorld;
  
  G4VPhysicalVolume* fPhysiVol;
  G4LogicalVolume*   fLogicVol;  
  G4Box*             fSolidVol;
  
  G4VPhysicalVolume* fPhysiGridVol;
  G4LogicalVolume*   fLogicGridVol;  
  G4Box*	     fSolidGridVol;

  G4VPhysicalVolume* fPhysiGridVol_Hole;
  G4LogicalVolume*   fLogicGridVol_Hole;  
  G4Box*	     fSolidGridVol_Hole;

  G4VPhysicalVolume* fPhysiControlVol_GridShadow;
  G4LogicalVolume*   fLogicControlVol_GridShadow;  
  G4Box*             fSolidControlVol_GridShadow;

  G4Material*        fDefaultMaterial;
  G4Material*        fGridMaterial;

  DetectorMessenger* fDetectorMessenger;
  
  //
  static G4ThreadLocal TabulatedField3D * fField;
  
};
#endif



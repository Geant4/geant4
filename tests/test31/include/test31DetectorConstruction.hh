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
#ifndef test31DetectorConstruction_h
#define test31DetectorConstruction_h 1

// -------------------------------------------------------------
//
//
// -------------------------------------------------------------
//      GEANT4 test31
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- test31DetectorConstruction -------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of test31 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Material.hh"
#include "G4UniformMagField.hh"
#include "globals.hh"
#include "G4ios.hh"

class test31DetectorMessenger;
class test31Histo;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class test31DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  test31DetectorConstruction();
  virtual ~test31DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();

  void SetAbsorberMaterial(const G4String& v);
  void SetGapMaterial(const G4String& v);
  void SetNumberOfAbsorbers(G4int);     
  void SetNumAbsorbersSaved(G4int);     
  void SetAbsorberThickness(G4double);     
  void SetAbsorberSizeXY   (G4double);            
  void SetWorldMaterial(const G4String& v);
  void SetWorldSizeZ   (G4double);
  void SetGap(G4double val);

  void SetFirstEventToDebug(G4int val) {nFirstEvtToDebug = val;};
  G4int GetFirstEventToDebug() const {return nFirstEvtToDebug;};
  void SetLastEventToDebug(G4int val) {nLastEvtToDebug = val;};
  G4int GetLastEventToDebug() const {return nLastEvtToDebug;};

  void SetVerbose(G4int val) {myVerbose = val;};
  G4int GetVerbose() const {return myVerbose;};

  void SetMagField(G4double,G4int);
     
  void UpdateGeometry();
  void PrintGeomParameters(); 
                                     
private:

  void DefineMaterials();
  void ComputeGeomParameters();

  void GeometryIsChanged();     

  G4String           nameMatWorld;
  G4String           nameMatAbsorber;
  G4String           nameMatGap;

  G4Material*        AbsorberMaterial;
  G4Material*        GapMaterial;
  G4Material*        WorldMaterial;
 
  G4double           WorldSizeZ;     
  G4double           SizeXY;
  G4double           AbsorberThickness;  // 
  G4double           gap;

  G4int              NumberOfAbsorbers;
  G4int              nSaved;

  G4Box*             solidWorld;    //pointer to the solid World 
  G4LogicalVolume*   logicWorld;    //pointer to the logical World
  G4VPhysicalVolume* physWorld;     //pointer to the physical World

  G4LogicalVolume*   logicGap;      //pointer to the logical Gap

  G4Box*             solidAbs;      //pointer to the solid Absorber
  G4LogicalVolume*   logicAbs;      //pointer to the logical Absorber
  G4VPhysicalVolume* physAbs;       //pointer to the physical Absorber
     
  G4UniformMagField* magField;      //pointer to the magnetic field
     
  test31DetectorMessenger* detectorMessenger;   
  test31Histo*       histo;

  G4int              myVerbose;      
  G4int              nFirstEvtToDebug;
  G4int              nLastEvtToDebug;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif


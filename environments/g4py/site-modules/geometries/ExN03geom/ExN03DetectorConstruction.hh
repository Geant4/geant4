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
// $Id: ExN03DetectorConstruction.hh,v 1.4 2006/06/29 15:29:18 gunter Exp $
// $Name: geant4-08-01 $
// ====================================================================
//   ExN03DetectorConstruction.hh
//
//                                         2005 Q
// ====================================================================
#ifndef EXN03_DETECTOR_CONSTRUCTION_H
#define EXN03_DETECTOR_CONSTRUCTION_H

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

// ====================================================================
//
// class definition
//
// ====================================================================
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;

class ExN03DetectorConstruction : public G4VUserDetectorConstruction {
private:
  G4Material*        AbsorberMaterial;
  G4double           AbsorberThickness;
  
  G4Material*        GapMaterial;
  G4double           GapThickness;
  
  G4int              NbOfLayers;
  G4double           LayerThickness;
  
  G4double           CalorSizeYZ;
  G4double           CalorThickness;
  
  G4Material*        defaultMaterial;
  G4double           WorldSizeYZ;
  G4double           WorldSizeX;
  
  G4Box*             solidWorld;    //pointer to the solid World 
  G4LogicalVolume*   logicWorld;    //pointer to the logical World
  G4VPhysicalVolume* physiWorld;    //pointer to the physical World
  
  G4Box*             solidCalor;    //pointer to the solid Calor 
  G4LogicalVolume*   logicCalor;    //pointer to the logical Calor
  G4VPhysicalVolume* physiCalor;    //pointer to the physical Calor
  
  G4Box*             solidLayer;    //pointer to the solid Layer 
  G4LogicalVolume*   logicLayer;    //pointer to the logical Layer
  G4VPhysicalVolume* physiLayer;    //pointer to the physical Layer
  
  G4Box*             solidAbsorber; //pointer to the solid Absorber
  G4LogicalVolume*   logicAbsorber; //pointer to the logical Absorber
  G4VPhysicalVolume* physiAbsorber; //pointer to the physical Absorber
  
  G4Box*             solidGap;      //pointer to the solid Gap
  G4LogicalVolume*   logicGap;      //pointer to the logical Gap
  G4VPhysicalVolume* physiGap;      //pointer to the physical Gap
  
  G4UniformMagField* magField;      //pointer to the magnetic field
  
  
  void DefineMaterials();
  void ComputeCalorParameters();
  G4VPhysicalVolume* ConstructCalorimeter();     

public:
  ExN03DetectorConstruction();
  ~ExN03DetectorConstruction();

  // set/get mehtods...
  void SetAbsorberMaterial (G4String);     
  void SetAbsorberThickness(G4double);     
  void SetGapMaterial (G4String);     
  void SetGapThickness(G4double);
  void SetCalorSizeYZ(G4double);          
  void SetNbOfLayers (G4int);   
  void SetMagField(G4double);
  
  G4double GetWorldSizeX()           { return WorldSizeX; }
  G4double GetWorldSizeYZ()          { return WorldSizeYZ; }
  G4double GetCalorThickness()       { return CalorThickness; } 
  G4double GetCalorSizeYZ()          { return CalorSizeYZ; }
  G4int GetNbOfLayers()              { return NbOfLayers; } 
  G4Material* GetAbsorberMaterial()  { return AbsorberMaterial; }
  G4double GetAbsorberThickness()    { return AbsorberThickness; } 
  
  G4Material* GetGapMaterial()       { return GapMaterial; }
  G4double GetGapThickness()         { return GapThickness; }
  
  const G4VPhysicalVolume* GetphysiWorld() { return physiWorld; }
  const G4VPhysicalVolume* GetAbsorber()   { return physiAbsorber; }
  const G4VPhysicalVolume* GetGap()        { return physiGap; }

  // operations...
  virtual G4VPhysicalVolume* Construct();
  void UpdateGeometry();
  void PrintCalorParameters();   

};


// ====================================================================
//   inline functions
// ====================================================================

inline void ExN03DetectorConstruction::ComputeCalorParameters() 
{
  // Compute derived parameters of the calorimeter
  LayerThickness = AbsorberThickness + GapThickness;
  CalorThickness = NbOfLayers*LayerThickness;
  
  WorldSizeX = 1.2*CalorThickness; WorldSizeYZ = 1.2*CalorSizeYZ;
}


#endif


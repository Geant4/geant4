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
// $Id: Tst50DetectorConstruction.hh,v 1.11 2003-05-28 08:10:10 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// author: Susanna Guatelli (guatelli@ge.infn.it) 
//
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
//
// -------------------------------------------------------------------

#ifndef Tst50DetectorConstruction_h
#define Tst50DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UserLimits;
class Tst50TrackerSD;
class Tst50DetectorMessenger;

class Tst50DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  Tst50DetectorConstruction();
  ~Tst50DetectorConstruction();   
  G4VPhysicalVolume* Construct();

private:
  void DefineMaterials();   
  G4VPhysicalVolume* ConstructWorld();     
     
public:
  void PrintParameters(); 
  void SetTargetMaterial (G4String); 
  G4double GetDensity(); //returns the slab absorber material density
  G4String  GetMaterialName(); //returns the absorber material name 
  void SetTargetThickness(G4double);    
  G4double  GetTargetThickness();    
  void SetTargetX(G4double); 
  void SetTargetY(G4double);  
  void UpdateGeometry(); 
  void SetMaxStepInTarget(G4double value); 
  void UseUserLimits(G4bool value); 
  void SetUserLimits(G4bool);
  G4Material* GetTargetMaterial()  {return targetMaterial;}; 
  //returns the absorber material

private: 
  G4bool  isRegisteredUserLimits;
  
  // available materials ...
  G4Material* hydrogen;
  G4Material* beryllium;
  G4Material* graphite; 
  G4Material* magnesium;
  G4Material* aluminium;
  G4Material* silicon;
  G4Material* liquidArgon;  
  G4Material* iron;   
  G4Material* gallium;
  G4Material* germanium;
  G4Material* silver;
  G4Material* cesium; 
  G4Material* gold; 
  G4Material* lead;
  G4Material* water; 
  G4Material* quartz; 
  G4Material* air; 
  G4Material* vacuum;
  
  G4Material*        targetMaterial;
  G4Material*        defaultMaterial;//World absorber material: vacuum
  
  G4Box*             solidWorld; 
  G4LogicalVolume*   logicWorld; 
  G4VPhysicalVolume* physiWorld; 
  G4Box*             solidTarget;
  G4LogicalVolume*   logicTarget;
  G4VPhysicalVolume* physiTarget;
  Tst50TrackerSD* targetSD; 
  
  G4double           targetThickness;
  G4double targetX;
  G4double targetY;
  G4UserLimits*    theUserLimitsForTarget;    
  G4bool           fUseUserLimits;
  G4double         theMaxStepInTarget;
  Tst50DetectorMessenger* messenger;
};
#endif





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
// $Id: Tst52DetectorConstruction.hh,v 1.1.2.1 2007-12-10 16:33:38 gunter Exp $
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

#ifndef Tst52DetectorConstruction_h
#define Tst52DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4NistManager.hh"
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UserLimits;
class Tst52TrackerSD;
class Tst52DetectorMessenger;
class Tst52PhantomROGeometry;

class Tst52DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  Tst52DetectorConstruction();
  ~Tst52DetectorConstruction();   
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
  void SetVoxelNumber(G4int);

private: 
  G4bool  isRegisteredUserLimits;
  // G4Material* water_nist;
  // available materials ...
  G4Material* hydrogen;
  G4Material* beryllium;
  G4Material* graphite; 
  G4Material* magnesium;
  G4Material* aluminium;
  G4Material* silicon;
  G4Material* liquidArgon;  
  G4Material* titanium;
  G4Material* iron; 
  G4Material* cobalt;
  G4Material* nickel;
  G4Material* indium;
  G4Material* tin; 
  G4Material* copper; 
  G4Material* zinc;  
  G4Material* gallium;
  G4Material* germanium;
  G4Material* zirconium;
  G4Material* molybdenium;
  G4Material* silver;
  G4Material* cadmium;
  G4Material* cesium; 
  G4Material* samarium;
  G4Material* ytterbium; 
  G4Material* tantalum;
  G4Material* tungsten;
  G4Material* gold; 
  G4Material* lead;
  G4Material* uranium;
  G4Material* water; 
  G4Material* quartz; 
  G4Material* ossigeno;
  G4Material* air; 
  G4Material* vacuum;
   G4Material* nytrogen;
  G4Material*        targetMaterial;
  G4Material*        defaultMaterial;//World absorber material: vacuum
  
  G4Box*             solidWorld; 
  G4LogicalVolume*   logicWorld; 
  G4VPhysicalVolume* physiWorld; 
  G4Box*             solidTarget;
  G4LogicalVolume*   logicTarget;
  G4VPhysicalVolume* physiTarget;
  Tst52TrackerSD* targetSD; 
  
  G4double           targetThickness;
  G4double targetX;
  G4double targetY;
  G4UserLimits*    theUserLimitsForTarget;    
  G4bool           fUseUserLimits;
  G4double         theMaxStepInTarget;
  Tst52DetectorMessenger* messenger; 
  Tst52PhantomROGeometry* targetROGeometry;
  // G4NistManager*     nistMan;
  G4int numberOfVoxelsAlongZ;
};
#endif





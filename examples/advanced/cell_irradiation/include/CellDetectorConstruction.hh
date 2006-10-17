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
//    **************************************
//    *                                    *
//    *    CellDetectorConstruction.hh     *
//    *                                    *
//    **************************************
//
//
// author: Susanna Guatelli (guatelli@ge.infn.it) 
//	   Barbara Mascialino (Barbara.Mascialino@ge.infn.it)
//
// History:
// -----------
// 20 September 2006   S. Guatelli, B. Mascialino   1st implementation
//
// -------------------------------------------------------------------

#ifndef CellDetectorConstruction_h
#define CellDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UserLimits;
class CellTrackerSD;
class CellDetectorMessenger;
class CellPhantomROGeometry;
class CellDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  CellDetectorConstruction();
  ~CellDetectorConstruction();   
  G4VPhysicalVolume* Construct();

private:
  // Definition of the materials 
  void DefineMaterials();   

  // Definition of the experimental set-up
  G4VPhysicalVolume* ConstructWorld();     
     
public:
 
  // Print the parameters of the target 
  void PrintParameters(); 
 
  // This method allows to change the material of 
  // target interactively
  void SetTargetMaterial (G4String); 

  G4double GetDensity(); // Returns the target material density
  G4String  GetMaterialName(); //returns the target material name 

  // This method allows to change the Z size
  // of the target interactively
  void SetTargetThickness(G4double);    
  G4double  GetTargetThickness(); // Returns the Z size of the target  

  // This method allows to change the X size
  // of the target interactively
  void SetTargetX(G4double); 

  // This method allows to change the X size
  // of the target interactively
  void SetTargetY(G4double);  

  void UpdateGeometry(); 
  
  // This method allows to change interactively 
  // the value of the maximum step in the target
  void SetMaxStepInTarget(G4double value); 

  void UseUserLimits(G4bool value); 
  void SetUserLimits(G4bool);

  G4Material* GetTargetMaterial()  {return targetMaterial;}; 
  // Returns the absorber material

  G4double GetTargetMass();


private: 
  G4bool  isRegisteredUserLimits;
  
  // Available materials ...
  G4Material* hydrogen;
  G4Material* beryllium;
  G4Material* graphite; 
  G4Material* magnesium;
  G4Material* aluminium; 
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
  CellTrackerSD* targetSD; 
  
  G4double           targetThickness;
  G4double targetX;
  G4double targetY;
  G4UserLimits*    theUserLimitsForTarget;    
  G4bool           fUseUserLimits;
  G4double         theMaxStepInTarget;
  CellDetectorMessenger* messenger;
  CellPhantomROGeometry* targetROGeometry;
};
#endif





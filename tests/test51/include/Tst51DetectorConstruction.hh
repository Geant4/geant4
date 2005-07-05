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
// $Id: Tst51DetectorConstruction.hh,v 1.1 2005-07-05 11:05:54 guatelli Exp $
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

#ifndef Tst51DetectorConstruction_h
#define Tst51DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Tubs;
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class Tst51DetectorMessenger;

class Tst51DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  Tst51DetectorConstruction();
  ~Tst51DetectorConstruction();   
  G4VPhysicalVolume* Construct();

private:
  void DefineMaterials();   
  G4VPhysicalVolume* ConstructWorld();     
     
public:
  void PrintParameters(); 
  void SetTargetMaterial (G4String); 
  void SetTargetThickness(G4double);    
  G4double  GetTargetThickness();    
  void UpdateGeometry(); 

  G4Material* GetTargetMaterial()  {return targetMaterial;}; 
  //returns the material of the target

private: 
  
  // available materials ...
  G4Material* hydrogen;
  G4Material* beryllium;
  G4Material* aluminium;
  G4Material* silicon;
  G4Material* iron; 
  G4Material* germanium;
  G4Material* gold; 
  G4Material* lead;
  G4Material* ossigeno;
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
  
  
  G4double targetThickness;
  G4double targetY;
  G4double targetX;
  Tst51DetectorMessenger* messenger;
};
#endif





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
// $Id: Mars01DetectorConstruction.hh,v 1.1 2001-12-13 14:58:41 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Mars01DetectorConstruction_h
#define Mars01DetectorConstruction_h 1

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class Mars01DetectorMessenger;

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class Mars01DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  Mars01DetectorConstruction();
  ~Mars01DetectorConstruction();
  
public:
  G4VPhysicalVolume* Construct();
  void SelectMaterial(G4String val);
  
  static const G4String& GetMaterialName(G4int idx);
  
private:
  enum { N_MaterialType = 3}; 
  enum { N_Regions      = 10}; 
  
  void SelectMaterialPointer();
  
  G4LogicalVolume*   simpleBoxLog[N_Regions];
  //  G4Material* materials[N_MaterialType][N_Regions];
  G4Material* materials[N_MaterialType];
  G4Material* selectedMaterial[N_Regions];
  
  static const    G4String materialName[N_MaterialType];
  static G4int    materialIndex;
  
  Mars01DetectorMessenger * detectorMessenger;    
  
};

#endif


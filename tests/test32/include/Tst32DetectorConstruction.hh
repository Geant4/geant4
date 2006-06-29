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
// $Id: Tst32DetectorConstruction.hh,v 1.2 2006-06-29 21:58:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst32DetectorConstruction_h
#define Tst32DetectorConstruction_h 1

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class Tst32DetectorMessenger;

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class Tst32DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  Tst32DetectorConstruction();
  ~Tst32DetectorConstruction();
  
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
  
  Tst32DetectorMessenger * detectorMessenger;    
  
};

#endif


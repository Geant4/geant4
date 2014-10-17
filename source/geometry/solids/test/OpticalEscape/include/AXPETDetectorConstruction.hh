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
// $Id$
// ------------------------------------------------------------
// Geant4 class header file
//
// This class is a class derived from G4VUserDetectorConstruction
// for constructing all particles and processes.
//
// 03/09/2008, by T.Nikitina
// ------------------------------------------------------------

#ifndef AXPETDetectorConstruction_h
#define AXPETDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class G4VSolid; 
class G4Box;
class G4Tubs;
class G4Sphere;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class AXPETMaterial;
class AXPETDetectorMessenger;

	
class AXPETDetectorConstruction: public G4VUserDetectorConstruction
{
 public:
 AXPETDetectorConstruction();
 ~AXPETDetectorConstruction();

 public:
 G4VPhysicalVolume* Construct();

 AXPETMaterial* pttoMaterial;
  // Solids 
  G4VSolid * GetSolid()  { return aVolume ; }
  void  SwitchDetector();
  G4VSolid* SelectDetector (const G4String& val);
  inline G4double GetExtend(){return lysoR;}
  inline void SetDetectorName(G4String Value){fval=Value;}
  // Rotation 
  inline void SetRotationInX(G4double value){xRot=value;}
  inline void SetRotationInY(G4double value){yRot=value;}
  inline void SetRotationInZ(G4double value){zRot=value;}
  inline G4double GetRotationInX(){return xRot;}
  inline G4double GetRotationInY(){return yRot;}
  inline G4double GetRotationInZ(){return zRot;}
  // Abortion of Run(for large statistics)
  inline void SetAbortAction(G4bool bval){fAbort=bval;}
  inline G4bool GetAbortAction(){return fAbort;}
 private: 
 G4Material* air;
 G4Material* lyso;
 G4Material* vacuum;
 
 G4Box*		       WorldVolume;
 G4LogicalVolume*      LogWorldVolume;
 G4VPhysicalVolume*    PhysWorldVolume;

 G4VSolid*	       LYSOVolume;
 G4LogicalVolume*      LogLYSOVolume;
 G4VPhysicalVolume*    PhysLYSOVolume;


 
 public:
 G4double worldXsize;
 G4double worldYsize;
 G4double worldZsize;

 G4double lysoR;
 G4double lysoDz;
 G4double lysoIz;
 
 private:
  
  AXPETDetectorMessenger* detectorMessenger;  // pointer to the Messenger
  
  G4String fval ;  
  G4VSolid* aVolume;
  G4double xRot;
  G4double yRot;
  G4double zRot;
  G4double fAbort;
};  

#endif

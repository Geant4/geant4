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
// $Id: BrachyDetectorConstruction.hh,v 1.3 2006/06/29 17:30:50 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//  
//    ****************************************
//    *                                      *
//    *    BrachyDetectorConstruction.hh     *
//    *                                      *
//    ****************************************

#ifndef BrachyDetectorConstruction_H
#define BrachyDetectorConstruction_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"

class BrachyPhantomSD;
class BrachyDetectorMessenger;
class G4LogicalVolume;
class G4Material;
class G4Tubs;
class G4Box;
class G4Sphere;
class G4Tubs;
class G4Colour;
class G4VPhysicalVolume;
class BrachyPhantomSD;
class BrachyPhantomROGeometry;
class G4VPhysicalVolume;
class BrachyMaterial;
class BrachyFactory;
class BrachyVoxelParameterisation;

class BrachyDetectorConstruction : public G4VUserDetectorConstruction
{
public:

  BrachyDetectorConstruction(G4String&);
  ~BrachyDetectorConstruction();

  G4VPhysicalVolume*   Construct();  
  void SwitchBrachytherapicSeed(); 
  void SelectBrachytherapicSeed(G4String val);
  void ConstructPhantom(); 
  void ConstructSensitiveDetector();
  void PrintDetectorParameters(); 
  void SetPhantomMaterial(G4String); 

  const G4double VoxelWidth_X(){return phantomDimensionX/numberOfVoxelsAlongX;}
  const G4double VoxelWidth_Z(){return phantomDimensionZ/numberOfVoxelsAlongZ;}
  const G4int   GetNumVoxelX()  {return  numberOfVoxelsAlongX;}
  const G4int   GetNumVoxelZ()  {return numberOfVoxelsAlongZ;}
  const G4double GetDimX()      {return phantomDimensionX;}
  const G4double GetBoxDim_Z()  {return  phantomDimensionZ;}

  void ComputeDimVoxel() {dimVoxel = phantomDimensionX/numberOfVoxelsAlongX;}

private:
  
  G4int detectorChoice; //Select brachytherapic seed
  BrachyPhantomSD* phantomSD;//pointer to sensitive detector
  BrachyPhantomROGeometry* phantomROGeometry;//pointer to ROGeometry
  BrachyFactory* factory;

  // World ...
  G4Box*             World;        //pointer to the solid World 
  G4LogicalVolume*   WorldLog;     //pointer to the logical World
  G4VPhysicalVolume* WorldPhys;    //pointer to the physical World

  // Phantom ... 
  G4Box*              Phantom;  //pointer to solid phantom
  G4LogicalVolume*    PhantomLog; //pointer to logic phantom
  G4VPhysicalVolume*  PhantomPhys; //pointer to physical phantom
  G4Material*         phantomAbsorberMaterial;
 
  G4double phantomDimensionX; //Phantom XDimension
  G4double phantomDimensionY; //Phantom YDimension
  G4double phantomDimensionZ; //Phantom ZDimension  
  G4int numberOfVoxelsAlongX; //Number of voxels along x axis
  G4int numberOfVoxelsAlongZ; //Number of voxels along z axis 
  G4double Worldx ; //World XDimension
  G4double Worldy ; //World YDimension
  G4double Worldz ; //World XDimension
  G4String sensitiveDetectorName; 
  BrachyDetectorMessenger* detectorMessenger; 
  BrachyMaterial* pMaterial; 
   
  G4double dimVoxel;   
};

#endif

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
// $Id: BrachyDetectorConstruction.hh,v 1.2 2004/05/25 08:36:17 guatelli Exp $
// GEANT4 tag $Name: geant4-06-02 $
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

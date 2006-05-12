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
// $Id: BrachyDetectorConstruction.hh,v 1.16 2006-05-12 13:23:48 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//  
//    ****************************************
//    *                                      *
//    *    BrachyDetectorConstruction.hh     *
//    *                                      *
//    ****************************************
// This class manages the geometry of the simulation experimental set-up
//

#ifndef BrachyDetectorConstruction_H
#define BrachyDetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"

class BrachyDetectorMessenger;
class G4LogicalVolume;
class G4Material;
class G4Box;
class G4Colour;
class G4VPhysicalVolume;
class BrachyPhantomSD;
class BrachyPhantomROGeometry;
class G4VPhysicalVolume;
class BrachyMaterial;
class BrachyFactory;

class BrachyDetectorConstruction : public G4VUserDetectorConstruction
{
public:

  BrachyDetectorConstruction(G4String&);
  ~BrachyDetectorConstruction();

  G4VPhysicalVolume*   Construct();  
  void SwitchBrachytherapicSeed(); //Change radiactive source through GUI
  void SelectBrachytherapicSeed(G4String val);
  void ConstructPhantom(); 
  void ConstructSensitiveDetector();
  void PrintDetectorParameters(); 
  void SetPhantomMaterial(G4String); 

  const G4double VoxelWidth_X(){return phantomSizeX/numberOfVoxelsAlongX;}
  const G4double VoxelWidth_Z(){return phantomSizeZ/numberOfVoxelsAlongZ;}
  const G4int    GetNumVoxelX(){return numberOfVoxelsAlongX;}
  const G4int    GetNumVoxelZ(){return numberOfVoxelsAlongZ;}
  const G4double GetDimX()     {return phantomSizeX;}
  const G4double GetBoxDim_Z() {return  phantomSizeZ;}

  void ComputeDimVoxel() {dimVoxel = phantomSizeX/numberOfVoxelsAlongX;}

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
 
  G4double phantomSizeX; //Phantom XSize
  G4double phantomSizeY; //Phantom YSize
  G4double phantomSizeZ; //Phantom ZSize  
  G4int numberOfVoxelsAlongX; //Number of voxels along x axis
  G4int numberOfVoxelsAlongY; //Number of voxels along y axis 
  G4int numberOfVoxelsAlongZ; //Number of voxels along z axis 
  G4double worldSizeX ; //World XSize
  G4double worldSizeY ; //World YSize
  G4double worldSizeZ ; //World XSize
  G4String sensitiveDetectorName; 
  BrachyDetectorMessenger* detectorMessenger; 
  BrachyMaterial* pMaterial; 
   
  G4double dimVoxel;   
};

#endif

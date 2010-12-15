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

  G4double VoxelWidth_X() const {return phantomSizeX/numberOfVoxelsAlongX;}
  G4double VoxelWidth_Z() const {return phantomSizeZ/numberOfVoxelsAlongZ;}
  G4int    GetNumVoxelX() const {return numberOfVoxelsAlongX;}
  G4int    GetNumVoxelZ() const {return numberOfVoxelsAlongZ;}
  G4double GetDimX()    const {return phantomSizeX;}
  G4double GetBoxDim_Z() const {return  phantomSizeZ;}

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

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

class G4VPhysicalVolume;



class G4VPhysicalVolume;
class BrachyDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  BrachyDetectorConstruction(G4String&);
  ~BrachyDetectorConstruction();


private:
  G4double m_BoxDimX;
  G4double m_BoxDimY;
  G4double m_BoxDimZ;
  G4int NumVoxelX;
  G4int NumVoxelZ;
  G4double dimVoxel;
  G4String m_SDName;

public:
  void PrintDetectorParameters(); 
  void SetAbsorberMaterial(G4String);
  const   G4double VoxelWidth_X() {return m_BoxDimX/NumVoxelX;};//num voxel
  const   G4double VoxelWidth_Z() {return m_BoxDimZ/NumVoxelZ;};
  const   G4int   GetNumVoxelX(){return NumVoxelX;};
  const   G4int   GetNumVoxelZ(){return NumVoxelZ;};
  const   G4double GetDimX(){return m_BoxDimX;};
       
  const   G4double GetBoxDim_Z() {return  m_BoxDimZ;};
  
 
 
private:

  BrachyDetectorMessenger* detectorMessenger; 

  G4Box*             ExpHall;    //pointer to the solid World 
  G4LogicalVolume*   ExpHallLog;    //pointer to the logical World
  G4VPhysicalVolume* ExpHallPhys;    //pointer to the physical World


  G4Box*              Phantom;
  G4LogicalVolume*    PhantomLog;
  G4VPhysicalVolume*   PhantomPhys;


  G4Tubs* Capsule ;
  G4LogicalVolume*  CapsuleLog;    //pointer to the logical World
  G4VPhysicalVolume* CapsulePhys;
 
  G4Sphere* CapsuleTip;
  G4LogicalVolume* CapsuleTipLog;
  G4VPhysicalVolume* CapsuleTipPhys;
   
  G4Tubs* IridiumCore;
  G4LogicalVolume* IridiumCoreLog;
  G4VPhysicalVolume* IridiumCorePhys;

  G4Material*              air;
  //G4Material*              water;
  G4Material*              CapsuleMat;
  G4Material*              IridiumMat;
  G4Material*            AbsorberMaterial;

  G4VPhysicalVolume* Construct();

private:
         
  void DefineMaterials();
  void ComputeDimVoxel();
   
public:
  G4VPhysicalVolume* ConstructDetector();  

 
};

inline void BrachyDetectorConstruction::ComputeDimVoxel()
{
  dimVoxel=m_BoxDimX/NumVoxelX; 

}

#endif









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
class BrachyFactoryIr;
class BrachyFactoryI;
class BrachyFactoryLeipzig;
class BrachyVoxelParameterisation;
class BrachyDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    BrachyDetectorConstruction(G4String&);
    ~BrachyDetectorConstruction();


private:
G4int* pVoxel;
G4double m_BoxDimX;
G4double m_BoxDimY;
G4double m_BoxDimZ;
G4int NumVoxelX;
G4int NumVoxelZ;
G4double dimVoxel;
G4String m_SDName;
G4int detectorChoice;
G4double ExpHall_x ;
G4double ExpHall_y ;
G4double ExpHall_z ;


BrachyPhantomSD* pPhantomSD;
BrachyPhantomROGeometry*   pPhantomROGeometry;


BrachyFactoryLeipzig* factoryLeipzig;
BrachyFactoryI* factoryI;
BrachyFactoryIr* factoryIr;
BrachyMaterial* pMat;


public:
void PrintDetectorParameters(); 
void SetAbsorberMaterial(G4String);

const   G4double VoxelWidth_X() {return m_BoxDimX/NumVoxelX;};//num voxel
const   G4double VoxelWidth_Z() {return m_BoxDimZ/NumVoxelZ;};
const   G4int   GetNumVoxelX(){return NumVoxelX;};
const   G4int   GetNumVoxelZ(){return NumVoxelZ;};
const   G4double GetDimX(){return m_BoxDimX;};

const   G4double GetBoxDim_Z() {return  m_BoxDimZ;};

void SwitchDetector();

private:

BrachyDetectorMessenger* detectorMessenger; 

G4Box*             ExpHall;    //pointer to the solid World 
G4LogicalVolume*   ExpHallLog;    //pointer to the logical World
G4VPhysicalVolume* ExpHallPhys;    //pointer to the physical World

G4Box*             solidVoxel;  
G4LogicalVolume*   logicVoxel;  
G4VPhysicalVolume* physiVoxel;  


G4Box*              Phantom;
G4LogicalVolume*    PhantomLog;
G4VPhysicalVolume*   PhantomPhys;
G4Box* Voxel;
G4LogicalVolume*  VoxelLog;
G4VPhysicalVolume* VoxelPhys;


G4Material*            AbsorberMaterial;


void  ConstructPhantom();
void ConstructSensitiveDetector();

private:
G4VPhysicalVolume*   Construct();  

 void ComputeDimVoxel();


public:

 G4VPhysicalVolume* ConstructDetector();
void SelectDetector(G4String val);



};

inline void BrachyDetectorConstruction::ComputeDimVoxel()
{
dimVoxel=m_BoxDimX/NumVoxelX; 

}

#endif









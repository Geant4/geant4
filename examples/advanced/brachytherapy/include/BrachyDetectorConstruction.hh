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



class BrachyWaterBoxSD;
class G4LogicalVolume;
class G4Material;
class G4Box;
class G4UserLimits;
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
         G4double dimVoxel;
	 G4int NumVoxelX;
	 G4int NumVoxelZ;
     
         G4String m_SDName;

public:
  void PrintDetectorParameters(); 
  void SetMaterial(G4String);
  void SetDimension(G4double);
  void SetNumVoxel(G4int);
  void UpdateGeometry();
  G4double VoxelWidth_X() {return m_BoxDimX/NumVoxelX;};//num voxel
  G4double VoxelWidth_Z() {return m_BoxDimZ/NumVoxelZ;};
  G4int   GetNumVoxelX(){return NumVoxelX;};
  G4int   GetNumVoxelZ(){return NumVoxelZ;};
  G4double GetDimX(){return m_BoxDimX;};
	
 // methods for UserLimits in Water
  void      UseUserLimits(G4bool value); 
  void  SetMaxStepInWater(G4double value); 

 G4bool           fUseUserLimits;
 G4UserLimits*    theUserLimitsForWater; 
 G4double         theMaxStepInWater;
 
 private:
 G4LogicalVolume*          WaterBoxLog;
  G4Material*              MaterialBox;
  G4Material*              DefaultMaterial;
  G4Material*              Iodio;
  G4Material*               Gold;
  G4Material*              titanium;
  
public:
	G4VPhysicalVolume* Construct();

 private:
         
     void DefineMaterials();
     void ComputeDimVoxel();
   

     G4VPhysicalVolume* ConstructDetector();  

 
};

inline void BrachyDetectorConstruction::ComputeDimVoxel()
{
 dimVoxel=m_BoxDimX/NumVoxelX; 

}

#endif







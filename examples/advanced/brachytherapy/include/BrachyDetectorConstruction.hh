//    ****************************************
//    *                                      *
//    *    BrachyDetectorConstruction.hh     *
//    *                                      *
//    ****************************************


#ifndef BrachyDetectorConstruction_H
#define BrachyDetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class BrachyWaterBoxSD;

class BrachyDetectorConstruction : public G4VUserDetectorConstruction
{
 public:
	BrachyDetectorConstruction(G4String &SDName,G4int NumVoxelX,G4int NumVoxelZ);
	~BrachyDetectorConstruction();

 public:
	const G4double m_BoxDimX;
	const G4double m_BoxDimY;
	const G4double m_BoxDimZ;

	const G4int m_NumVoxelX;
	const G4int m_NumVoxelZ;

	G4String m_SDName;

 public:
	G4VPhysicalVolume* Construct();
};
#endif


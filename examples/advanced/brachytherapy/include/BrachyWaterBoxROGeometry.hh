//    ************************************
//    *                                  *
//    *    BrachyWaterBoxROGeometry.hh   *
//    *                                  *
//    ************************************


#ifndef BrachyWaterBoxROGeometry_h
#define BrachyWaterBoxROGeometry_h 1

#include "G4VReadOutGeometry.hh"

class BrachyWaterBoxROGeometry : public G4VReadOutGeometry
{
  public:
	BrachyWaterBoxROGeometry(G4String aString,G4double DetDimX,G4double DetDimZ,G4int NumVoxelX,G4int NumVoxelZ);
  	~BrachyWaterBoxROGeometry();

  public:
	const G4double m_DetDimX;
	const G4double m_DetDimZ;
	const G4int m_NumVoxelX;
	const G4int m_NumVoxelZ;

  private:
  	G4VPhysicalVolume* Build();
};

#endif

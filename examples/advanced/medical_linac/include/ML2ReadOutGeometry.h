#ifndef CML2ReadOutGeometryH
#define CML2ReadOutGeometryH

#include "G4VReadOutGeometry.hh"

class G4VPhysicalVolume;

class CML2ReadOutGeometry : public G4VReadOutGeometry
{
public:
	CML2ReadOutGeometry(const G4RotationMatrix *m, G4ThreeVector *v);
	~CML2ReadOutGeometry(void);
	void setBuildData(G4ThreeVector centre, G4ThreeVector halfSize, G4int NumberOfVoxelsAlongX, G4int NumberOfVoxelsAlongY, G4int NumberOfVoxelsAlongZ);
	G4VPhysicalVolume* Build();
private:
	G4VPhysicalVolume *ROPhyVol;
	G4ThreeVector centre, halfSize;
	G4int NumberOfVoxelsAlongX, NumberOfVoxelsAlongY, NumberOfVoxelsAlongZ;
};

#endif


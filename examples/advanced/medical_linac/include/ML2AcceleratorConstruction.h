#ifndef CML2AcceleratorConstructionH
#define CML2AcceleratorConstructionH

#include "ML2SinputData.h"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "ML2Acc1.h"

class CML2Acc1;
class CML2AcceleratorConstructionMessenger;
class CML2AcceleratorConstruction
{
public:
	CML2AcceleratorConstruction(void);
	~CML2AcceleratorConstruction(void);
	static CML2AcceleratorConstruction* GetInstance(void);
	void Construct(G4VPhysicalVolume *PVWorld);
	G4VPhysicalVolume *getPhysicalVolume(void){return this->PVAccWorld;};
	inline void setAcceleratorName(G4String val){this->AcceleratorName=val;};
	inline void setAcceleratorSpecficationsFileName(G4String val){this->AcceleratorSpecficationsFileName=val;};
	inline void setAcceleratorRotationX(G4double val){this->rotationX=val;};
	inline void setAcceleratorRotationY(G4double val){this->rotationY=val;};
	inline void setAcceleratorRotationZ(G4double val){this->rotationZ=val;};
private:
	void design(void);
	CML2AcceleratorConstructionMessenger *acceleratorConstructionMessenger;
	static CML2AcceleratorConstruction * instance;
	G4String AcceleratorName, AcceleratorSpecficationsFileName;

	G4VPhysicalVolume *PVAccWorld;
	G4RotationMatrix *rotation;
	G4double rotationX, rotationY, rotationZ;
	CML2Acc1 *accelerator1;
};


#endif

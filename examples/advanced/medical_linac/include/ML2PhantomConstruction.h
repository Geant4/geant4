#ifndef CML2PhantomConstructionH
#define CML2PhantomConstructionH

#include "globals.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"

#include "ML2SinputData.h"
#include "ML2Ph_FullWater.h"
#include "ML2Ph_BoxInBox.h"

class CML2PhantomConstructionMessenger;
class CPh_HalfWaterAir;
class CPh_FullWater;
class CML2PhantomConstruction
{
public:
	CML2PhantomConstruction(void);
	~CML2PhantomConstruction(void);
	static CML2PhantomConstruction* GetInstance(void);
	void Construct(G4VPhysicalVolume *PVWorld, G4int saving_in_ROG_Voxels_every_events, G4int seed, G4String ROGOutFile, G4bool bSaveROG);
	inline G4int getTotalNumberOfEvents()
	{
		if (this->phantomName="fullWater")
		{return this->Ph_fullWater->getTotalNumberOfEvents();}
		else if (this->phantomName="BoxInBox")
		{return this->Ph_BoxInBox->getTotalNumberOfEvents();}
		return 0;
	};
	inline void setPhantomName(G4String val){this->phantomName=val;};
	inline void setPhantom_nVoxelsX(G4int val){this->nVoxelsX=val;};
	inline void setPhantom_nVoxelsY(G4int val){this->nVoxelsY=val;};
	inline void setPhantom_nVoxelsZ(G4int val){this->nVoxelsZ=val;};
	inline void setPhantomSpecficationsFileName(G4String val){this->PhantomSpecficationsFileName=val;};
	inline void setPhantomRotationX(G4double val){this->rotationX=val;};
	inline void setPhantomRotationY(G4double val){this->rotationY=val;};
	inline void setPhantomRotationZ(G4double val){this->rotationZ=val;};
private:
	void design(void);
	CML2PhantomConstructionMessenger *phantomContstructionMessenger;
	static CML2PhantomConstruction * instance;
	G4int nVoxelsX, nVoxelsY, nVoxelsZ;
	G4String phantomName, PhantomSpecficationsFileName;

	G4RotationMatrix *rotation;
	G4VPhysicalVolume *PVPhmWorld;

	G4double rotationX, rotationY, rotationZ;
	CML2Ph_FullWater *Ph_fullWater;
	CML2Ph_BoxInBox *Ph_BoxInBox;
};
#endif


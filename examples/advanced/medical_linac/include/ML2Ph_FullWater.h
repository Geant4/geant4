#ifndef CML2Ph_FullWaterH
#define CML2Ph_FullWaterH

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"

#include "ML2SDWithParticle.h"
#include "ML2SDWithVoxels.h"
#include "ML2ReadOutGeometry.h"

#include "G4SDManager.hh"
#include "G4ProductionCuts.hh"


class CML2Ph_FullWater
{
public:
	CML2Ph_FullWater();
	~CML2Ph_FullWater(void);
	bool Construct(G4VPhysicalVolume *PVWorld, G4int saving_in_ROG_Voxels_every_events, G4int seed, G4String ROGOutFile, G4bool bSaveROG);
	inline G4int getTotalNumberOfEvents(){return this->sensDet->getTotalNumberOfEvents();};
private:
	G4VPhysicalVolume *PVWorld;
	CML2SDWithVoxels *sensDet;
};


#endif

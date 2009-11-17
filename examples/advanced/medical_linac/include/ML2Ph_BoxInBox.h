#ifndef CML2Ph_BoxInBoxH
#define CML2Ph_BoxInBoxH

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "ML2SDWithParticle.h"
#include "ML2SDWithVoxels.h"
#include "ML2ReadOutGeometry.h"
#include "G4BooleanSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4ProductionCuts.hh"
#include "G4UserLimits.hh"

class CML2Ph_BoxInBox
{
public:
	CML2Ph_BoxInBox();
	~CML2Ph_BoxInBox(void);
	bool Construct(G4VPhysicalVolume *PVWorld, G4int saving_in_ROG_Voxels_every_events, G4int seed, G4String ROGOutFile, G4bool bSaveROG);
	inline G4int getTotalNumberOfEvents(){return this->sensDet->getTotalNumberOfEvents();};
private:
	G4VPhysicalVolume *PVWorld;
	CML2SDWithVoxels *sensDet;

	G4ThreeVector centreBoxInside;
	G4double halfBoxInside_Thickness; 
};

#endif


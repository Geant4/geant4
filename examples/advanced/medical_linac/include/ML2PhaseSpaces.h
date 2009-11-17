#ifndef CML2PhaseSpacesH
#define CML2PhaseSpacesH

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"

#include "ML2SDWithParticle.h"
#include "ML2SDWithVoxels.h"
#include "G4SDManager.hh"
#include "ML2SinputData.h"

class CML2PhaseSpaces
{
public:
	CML2PhaseSpaces();
	~CML2PhaseSpaces(void);
	bool createPlane(G4VPhysicalVolume  *PVWorld, G4String name, G4ThreeVector centre, G4ThreeVector halfSize);
	bool createPlane(G4int idSD_Type, G4int max_N_particles_in_PhSp_File, G4int seed, G4int nMaxParticlesInRamPhaseSpace, G4VPhysicalVolume  *PVWorld, G4String name, G4String PhaseSpaceOutFile, G4bool bSavePhaseSpace, G4bool bStopAtPhaseSpace, G4ThreeVector centre, G4ThreeVector halfSize, SPrimaryParticle *primaryParticleData);
	G4int getCML2SensDetNParticle(){return this->sensDetParticle->getTotalNumberOfParticles();};
	inline CML2SDWithParticle* getCML2SensitiveDetectorParticle(){return this->sensDetParticle->getCML2SensitiveDetectorParticle();};
	inline bool getBContinueRun(){return this->sensDetParticle->getBContinueRun();};
	CML2SDWithParticle *sensDetParticle;
private:
	CML2SDWithVoxels *sensDetVoxelized;
	G4int nParticles;
};


#endif

#ifndef CCML2SDWithVoxelsH
#define CCML2SDWithVoxelsH

#include "G4VSensitiveDetector.hh"
#include "G4TouchableHistory.hh"
#include "G4VTouchable.hh"
#include "ML2SinputData.h"

class CML2ReadOutGeometryVoxels;


class CML2SDWithVoxels : public G4VSensitiveDetector
{
public:
	CML2SDWithVoxels(G4String name, G4int saving_in_ROG_Voxels_every_events, G4int seed, G4String ROGOutFile, G4bool bSaveROG, G4ThreeVector halfSize, G4int NumberOfVoxelsAlongX, G4int NumberOfVoxelsAlongY, G4int NumberOfVoxelsAlongZ);
	~CML2SDWithVoxels(void);
	G4bool ProcessHits(G4Step *aStep, G4TouchableHistory *ROHist);
	inline void Initialize(G4HCofThisEvent *){};
	inline void EndOfEvent(G4HCofThisEvent*){};
	G4int getTotalNumberOfEvents(){return this->nTotalEvents;};
	inline void setActive(G4bool bActive){this->bActive=bActive;};

private:
	void saveDataInVoxels();
	void saveHeaderDataInVoxels();

	G4ThreeVector halfSize;
	G4ThreeVector pos;
	G4double halfXVoxelDimensionX, halfXVoxelDimensionY, halfXVoxelDimensionZ;
	G4int NumberOfVoxelsAlongX, NumberOfVoxelsAlongY, NumberOfVoxelsAlongZ;
	Svoxel ***voxels;
	G4String fullOutFileData;
	G4int nTotalEvents, nParticle, nParticleValatile, saving_in_ROG_Voxels_every_events;
	G4bool bActive, bSaveROG;
	G4double voxelMass, density, voxelVolume;
	
};



#endif

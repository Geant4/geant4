#ifndef CML2SensitiveDetectorParticleH
#define CML2SensitiveDetectorParticleH

#include "G4VSensitiveDetector.hh"
#include "G4TouchableHistory.hh"
#include "G4VTouchable.hh"
#include "ML2SinputData.h"

class CML2ReadOutGeometryVoxels;

class CML2SDWithParticle : public G4VSensitiveDetector
{
public:
	CML2SDWithParticle();
	CML2SDWithParticle(G4int idType, G4int max_N_particles_in_PhSp_File, G4int seed, G4int nMaxParticlesInRamPhaseSpace, G4String name, G4String PhaseSpaceOutFile, G4bool bSavePhaseSpace, G4bool bStopAtVolatilePhaseSpace, SPrimaryParticle *primaryParticleData);
	~CML2SDWithParticle(void);
	G4bool ProcessHits(G4Step *aStep, G4TouchableHistory *ROHist);
//	inline void setbStopAtPhaseSpace(G4bool bStopAtPhaseSpace){this->bStopAtPhaseSpace=bStopAtPhaseSpace;};
	//inline void Initialize(G4HCofThisEvent *){};
	//inline void EndOfEvent(G4HCofThisEvent*){};
	G4int getTotalNumberOfParticles(){return this->nTotalParticles;};
	inline CML2SDWithParticle* getCML2SensitiveDetectorParticle(){return this;};
	inline Sparticle getParticle(int i){return this->particles[i];};
	inline void setActive(G4bool bActive){this->bActive=bActive;};
	inline bool getBContinueRun(){return this->bContinueRun;};
private:
	void saveDataParticles(G4int nParticle);
	void saveHeaderParticles();

	G4ThreeVector halfSize;
	G4ThreeVector pos;
	SPrimaryParticle *primaryParticleData;
	Sparticle *particles;
	G4String fullOutFileData;
	G4int nTotalParticles, nParticle;
	G4int idType, nMaxParticlesInRamPhaseSpace, max_N_particles_in_PhSp_File;
	G4bool bStopAtPhaseSpace, bSavePhaseSpace, bActive, bContinueRun;
};

#endif

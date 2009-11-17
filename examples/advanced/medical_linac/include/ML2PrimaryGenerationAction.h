#ifndef CML2PrimaryGenerationActionH
#define CML2PrimaryGenerationActionH


#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "G4Timer.hh"
#include "Randomize.hh" 

#include "G4ParticleDefinition.hh"
#include "ML2SinputData.h"
#include "ML2SDWithParticle.h"
#include "ML2SDWithVoxels.h"

class G4ParticleGun;
class G4ParticleDefinition;
class CML2PrimaryGenerationActionMessenger;

class CML2PrimaryGenerationAction : public G4VUserPrimaryGeneratorAction
{
public:
	CML2PrimaryGenerationAction(SPrimaryParticle *primaryParticleData);
	~CML2PrimaryGenerationAction(void);
	void design();
	void GeneratePrimaries(G4Event *anEvent);
	inline void setNIdenticalParticles(G4int val){this->nIdenticalParticles=val;};
	inline void setNLoopsPhSpParticles(G4int val){this->nLoopsPhSpParticles=val;};
	inline void setNMaxParticlesInRamPhaseSpace(G4int val){this->nMaxParticlesInRamPhaseSpace=val;};

	inline void setGunMeanEnergy(G4double val){this->GunMeanEnegy=val;};
	inline void setGunStdEnergy(G4double val){this->GunStdEnegy=val;};
	inline void setGunRadious(G4double val){this->GunRadious=val;};
	inline void setCalculatedPhaseSpaceFileIN(G4String val){this->calculatedPhaseSpaceFileIN=val;};
	inline void setSourceTypeName(G4String val)
	{
		this->sourceTypeName=val;
		if (this->sourceTypeName=="randomTarget")
		{
			this->idParticleSource=id_randomTarget;
		}
		else if (this->sourceTypeName=="phaseSpace")
		{
			this->idParticleSource=id_phaseSpace;
		}
	};
	
private:
	void setGunRandom();
	void setGunCalculatedPhaseSpace();
	void GenerateFromRandom();
	void GenerateFromCalculatedPhaseSpace();
	void fillParticlesContainer();
	bool itIsTheSameParticle(Sparticle *p1, Sparticle *p2);

	G4int nBeam, nIdenticalParticles, nLoopsPhSpParticles, idGunType, nMaxParticlesInRamPhaseSpace, idParticleSource;
	G4double GunMeanEnegy, GunStdEnegy, GunRadious;
	G4String calculatedPhaseSpaceFileIN;

	CML2PrimaryGenerationActionMessenger *PrimaryGenerationActionMessenger;


	G4ThreeVector dir, pos;
	G4double ek;

	G4Timer myTime;
	G4double sinTheta, cosTheta, phi;
	G4double ro, alfa;
	G4ParticleGun *particleGun;
	G4ParticleDefinition *gamma;
	G4ParticleDefinition *electron;
	G4ParticleDefinition *positron;
	G4int nEventsPerRun;
	SPrimaryParticle *primaryParticleData;
	Sparticle *particles, *particle;
	int nParticle, nPhSpParticles, nRandomParticles, idCurrentParticleSource;
	Sparticle *firstFileParticle, *lastLoadedParticle;
	G4String sourceTypeName;
};

#endif

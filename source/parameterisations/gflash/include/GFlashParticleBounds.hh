// Created by Joanna Weng 9.11.2004 

#ifndef GFlashParticleBounds_h
#define GFlashParticleBounds_h 

#include  "G4ParticleDefinition.hh"

class  GFlashParticleBounds
{
	public:
	GFlashParticleBounds();
	~GFlashParticleBounds();
	
	// methods to get/set ELE/Gamma max & min energy bounds
	G4double GetMinEneToParametrise(G4ParticleDefinition &particleType);
	G4double GetMaxEneToParametrise(G4ParticleDefinition &particleType); 
	G4double GetEneToKill(G4ParticleDefinition &particleType) ;
	
	void SetMinEneToParametrise(G4ParticleDefinition &particleType,G4double enemin);
	void SetMaxEneToParametrise(G4ParticleDefinition &particleType,G4double enemax);
	void SetEneToKill(G4ParticleDefinition &particleType,G4double enekill);
	private:
	
	// electron and positron
	G4double EMinEneToParametrise;
	G4double EMaxEneToParametrise;
	G4double EEneToKill;
};
#endif


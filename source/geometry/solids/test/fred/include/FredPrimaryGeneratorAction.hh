//
// FredPrimaryGeneratorAction
//
// Define Fred's generator
//

#ifndef FredPrimaryGeneratorAction_h
#define FredPrimaryGeneratorAction_h

#include "G4VUserPrimaryGeneratorAction.hh"

#include "FredMessenger.hh"

class SprayParticleGun;
class GridParticleGun;
class G4ParticleGun;
class G4Event;

class FredPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
	public:
	FredPrimaryGeneratorAction( FredMessenger *ourMessenger );
	~FredPrimaryGeneratorAction();
	
	void GeneratePrimaries( G4Event *anEvent );
	
	private:
	SprayParticleGun *sprayGun;
	GridParticleGun  *gridGun;
	G4ParticleGun	 *g4Gun;
	
	FredMessenger	*messenger;
};

#endif

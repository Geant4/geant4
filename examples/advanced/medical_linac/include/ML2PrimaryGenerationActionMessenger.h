#ifndef ML2PrimaryGenerationActionMessengerH
#define ML2PrimaryGenerationActionMessengerH


#include "globals.hh"
#include "G4UImessenger.hh"

class CML2PrimaryGenerationAction;
class G4UImessenger;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

class CML2PrimaryGenerationActionMessenger : public G4UImessenger 
{
public:
	CML2PrimaryGenerationActionMessenger(CML2PrimaryGenerationAction *PML2PrimaryGenerationAction);
	~CML2PrimaryGenerationActionMessenger(void);
	void SetNewValue(G4UIcommand* cmd, G4String newValue);

private:
	CML2PrimaryGenerationAction *pML2PrimaryGenerationAction;

	G4UIcmdWithAnInteger *nIdenticalParticles, *nLoopsPhSpParticles;
	G4UIcmdWithAnInteger *nMaxParticlesInRamPhaseSpace;
	G4UIcmdWithADoubleAndUnit *GunMeanEnegy, *GunStdEnegy, *GunRadious;
	G4UIcmdWithAString  *calculatedPhaseSpaceFileIN, *sourceTypeName;
	G4UIcmdWithADoubleAndUnit *rotationX,  *rotationY,  *rotationZ;
};

#endif


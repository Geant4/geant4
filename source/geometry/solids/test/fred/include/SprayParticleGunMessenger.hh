//
// SprayParticleGunMessenger.hh
//
// Declaration of SprayParticleGun's UI interface
//

#ifndef SprayParticleGunMessenger_HH
#define SprayParticleGunMessenger_HH

class SprayParticleGun;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;

#include "G4UImessenger.hh"
#include "globals.hh"

class SprayParticleGunMessenger : public G4UImessenger
{
	public:
	SprayParticleGunMessenger( SprayParticleGun *gun );
	~SprayParticleGunMessenger();
	
	void SetNewValue( G4UIcommand *command, G4String newValues );
	G4String GetCurrentValue( G4UIcommand *command );
	
	private:
	SprayParticleGun	*gun;
	
	G4UIdirectory			*gunDirectory;
	G4UIcmdWith3VectorAndUnit	*positionCmd;
	G4UIcmdWithAnInteger		*xSprayCmd;
	G4UIcmdWithAnInteger		*ySprayCmd;
	G4UIcmdWithAnInteger		*zSprayCmd;
};

#endif

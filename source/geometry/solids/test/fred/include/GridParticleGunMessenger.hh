//
// GridParticleMessenger.hh
//
// Declaration of GridParticleGun's messenger
//
#ifndef GridParticleGunMessenger_hh
#define GridParticleGunMessenger_hh

class GridParticleGun;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;

#include "G4UImessenger.hh"
#include "globals.hh"

class GridParticleGunMessenger : public G4UImessenger
{
	public:
	GridParticleGunMessenger( GridParticleGun *gun );
	~GridParticleGunMessenger();
	
	void SetNewValue( G4UIcommand *command, G4String newValues );
	G4String GetCurrentValue( G4UIcommand *command );
	
	private:
	GridParticleGun	*gun;
	
	G4UIdirectory				*gunDirectory;
	G4UIcmdWith3Vector			*directionCmd;
	G4UIcmdWith3VectorAndUnit	*originCmd;
	G4UIcmdWith3VectorAndUnit	*grid1Cmd;
	G4UIcmdWith3VectorAndUnit	*grid2Cmd;
	G4UIcmdWithAnInteger		*n1Cmd;
	G4UIcmdWithAnInteger		*n2Cmd;
};

#endif

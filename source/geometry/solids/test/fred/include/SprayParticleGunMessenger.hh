//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
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

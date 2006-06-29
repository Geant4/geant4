//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
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

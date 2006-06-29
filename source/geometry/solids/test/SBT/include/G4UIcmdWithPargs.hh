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
// G4UIcmdWithPargs.hh
//
// An extended type of G4UIcommand that accepts a run-time
// specified list of generic argument objects.
//

#ifndef G4UIcmdWithPArgs_hh
#define G4UIcmdWithPArgs_hh

#include "G4UIcommand.hh"

class G4UIcmdParg;

class G4UIcmdWithPargs : public G4UIcommand {
	public:
	G4UIcmdWithPargs( const char *path, G4UImessenger *messenger,
			  G4UIcmdParg **args, const int numParg );
	virtual ~G4UIcmdWithPargs();
	
	G4bool GetArguments( G4String argumentString );
	
	protected:
	G4UIcmdParg	**args;
	G4int		numArgs;
};

#endif

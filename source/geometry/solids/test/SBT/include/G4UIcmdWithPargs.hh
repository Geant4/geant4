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

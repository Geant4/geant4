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

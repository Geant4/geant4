//
// G4UIcmdParg.hh
//
// Specification of a pure virtual class of generic command
// arguments, to be used by G4UIcmdWithPargs.
//

#ifndef G4UIcmdParg_hh
#define G4UIcmdParg_hh

#include "globals.hh"

class istream;

class G4UIcmdParg {
	public: 
	G4UIcmdParg( const G4String &theName ) { name = theName; }
	virtual ~G4UIcmdParg() {;}
	
	virtual G4String ConvertToString() = 0;
	
	virtual char GetTypeCode() const = 0;
	inline G4String GetName() const { return name; }
	
	virtual istream &FetchValue( istream &ios ) = 0;
	
	protected:
	G4String name;
};


#endif

//
// G4UIcmdWithPargs.cc
//
// Implementation of a G4UIcommand with generic arguments
//

#include "G4UIcmdWithPargs.hh"
#include "G4UIcmdParg.hh"

//
// Constructor
//
G4UIcmdWithPargs::G4UIcmdWithPargs( const char *path, G4UImessenger *theMessenger,
			  	    G4UIcmdParg **theArgs, const int theNumParg )
		: G4UIcommand( path, theMessenger )
{
	args = theArgs;
	numArgs = theNumParg;
	
	G4UIcmdParg **thisArg;
	G4int index;
	for( index=0, thisArg=args; index < numArgs; index++, thisArg++ ) {
		//
		// Boy, sure looks like a memory leak, but it isn't.
		// See G4UIcommand::~G4UIcommand
		//
		G4UIparameter *param = new G4UIparameter( (*thisArg)->GetTypeCode() );
		param->SetParameterName( (*thisArg)->GetName() );
		param->SetOmittable(false);
		SetParameter( param );
	}
}


//
// Destructor
//
G4UIcmdWithPargs::~G4UIcmdWithPargs() {;}


//
// GetArguments
//
// Unlike the standard G4UIcommand, the values are returned in the
// in the array of arguments given to the constructor
//
G4bool G4UIcmdWithPargs::GetArguments( G4String argumentString )
{
	const char *buff = argumentString;
	istrstream is( (char *)buff );
	
	G4UIcmdParg **thisArg = args;
	while( thisArg < args + numArgs ) {
		(*thisArg++)->FetchValue(is);
		if (is.fail()) return false;
	}
	
	return true;
}

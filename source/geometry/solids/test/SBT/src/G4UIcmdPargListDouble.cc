//
// G4UIcmdPargListDouble.cc
//
// Implementation of an argument consisting of
// a list of double values for G4UIcmdWithParg
//

#include "G4UIcmdPargListDouble.hh"

#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif

//
// Constructor
//
G4UIcmdPargListDouble::G4UIcmdPargListDouble( const G4String &theName, G4int theMaxItem, 
					      const G4double theUnits )
			: G4UIcmdPargList( theName, theMaxItem )
{
	units = theUnits;
	//
	// Allocate memory for maximum number of items
	//
	storage = new G4double [theMaxItem];
}


//
// Destructor
//
G4UIcmdPargListDouble::~G4UIcmdPargListDouble()
{
	delete [] storage;
}


//
// FetchItem
//
G4bool G4UIcmdPargListDouble::FetchItem( const char *string, const G4int item )
{
	istrstream is( (char *)string );
	is >> storage[item];
	if (is.fail()) {
		G4cerr << "Error interpreting '" << string << "' as a double value" << endl;
		return false;
	}
	storage[item] *= units;
	
	return true;
}


//
// WriteItem
//
ostream &G4UIcmdPargListDouble::WriteItem( ostream &ios, const G4int item )
{
	ios << storage[item];
	return ios;
}

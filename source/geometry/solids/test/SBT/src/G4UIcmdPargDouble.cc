//
// G4UIcmdPargDouble.cc
//
// Implementation of a double argument for G4UIcmdWithPargs
//

#include "G4UIcmdPargDouble.hh"

#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif

//
// Constructor
//
G4UIcmdPargDouble::G4UIcmdPargDouble( const G4String &theName, 
				      const G4double def, const G4double theUnits )
			: G4UIcmdParg( theName ) 
{
	value = def;
	units = theUnits;
}


//
// Fetch argument value from input stream
//
istream &G4UIcmdPargDouble::FetchValue( istream &istr )
{
	G4double newValue;
	istr >> newValue;
	
	if (istr) value = newValue*units;
	
	return istr;
}


//
// ConvertToString
//
G4String G4UIcmdPargDouble::ConvertToString()
{
	char buff[20];
	ostrstream os(buff,20);
	
	os << value << '\0';
	
	G4String answer(buff);
	return answer;
}

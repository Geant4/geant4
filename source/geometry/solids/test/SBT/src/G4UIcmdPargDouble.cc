//
// G4UIcmdPargDouble.cc
//
// Implementation of a double argument for G4UIcmdWithPargs
//

#include "G4UIcmdPargDouble.hh"

#include "g4std/strstream"

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
G4std::istream &G4UIcmdPargDouble::FetchValue( G4std::istream &istr )
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
	G4std::ostrstream os(buff,20);
	
	os << value << '\0';
	
	G4String answer(buff);
	return answer;
}

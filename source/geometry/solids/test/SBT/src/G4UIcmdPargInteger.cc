//
// G4UIcmdPargInteger.cc
//
// Implementation of an integer argument for G4UIcmdWithPargs
//

#include "G4UIcmdPargInteger.hh"

#include "g4std/strstream"

//
// Fetch argument value from input stream
//
G4std::istream &G4UIcmdPargInteger::FetchValue( G4std::istream &istr )
{
	G4int newValue;
	istr >> newValue;
	
	if (istr) value = newValue;
	
	return istr;
}


//
// ConvertToString
//
G4String G4UIcmdPargInteger::ConvertToString()
{
	char buff[20];
	G4std::ostrstream os(buff,20);
	
	os << value << '\0';
	
	G4String answer(buff);
	return answer;
}

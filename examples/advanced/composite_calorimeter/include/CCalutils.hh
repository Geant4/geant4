///////////////////////////////////////////////////////////////////////////////
// File: CCalutils.hh
// Description: General utilities.
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalutils_hh
#define CCalutils_hh 1

#include "g4std/iostream"
#include "g4std/fstream"
#include "globals.hh"


// "number " + i = "number i"
G4String operator+(const G4String&, const int);
G4String operator+(const G4String&, const double);

//readName(G4std::ifstream& is, G4String& name) reads a name into G4String between 
//                                             quotes and skips lines begining 
//                                             with #. and if found *ENDDO returns
G4std::ifstream& readName(G4std::ifstream&, G4String&);

//findDO(G4std::ifstream& is, const G4String& str) reads until a *DO str is found
G4std::ifstream& findDO(G4std::ifstream&, const G4String&);

//tabs
G4std::ostream& tab(G4std::ostream&);

//ignores character until end of line
G4std::istream& jump(G4std::istream&);

//Opens the geometry file, either locally (if it exists) or "remotely"
bool openGeomFile(G4std::ifstream& is, const G4String& pathname, const G4String& filename);

//extracts value from a number*dimension (e.g. 10*mm)
G4double getDoubleValue( G4String paramString );

#endif

///////////////////////////////////////////////////////////////////////////////
// File: utils.hh
// Date: 03/98 I. Gonzalez
// Modifications: 29/03/98 I.G.
//                27/03/00 S.B. In OSCAR
// Description: General utilities.
///////////////////////////////////////////////////////////////////////////////
#ifndef utils_hh
#define utils_hh 1

#include <iostream>
#include <fstream>
#include "globals.hh"


// "number " + i = "number i"
G4String operator+(const G4String&, const int);
G4String operator+(const G4String&, const double);

//readName(istream& is, G4String& name) reads a name into G4String between 
//                                      quotes and skips lines begining 
//                                      with #. and if found *ENDDO returns
istream& readName(istream&, G4String&);

//findDO(istream& is, const G4String& str) reads until a *DO str is found
istream& findDO(istream&, const G4String&);

//tabs
ostream& tab(ostream&);

//ignores character until end of line
istream& jump(istream&);

//Opens the geometry file, either locally (if it exists) or "remotely"
bool openGeomFile(ifstream& is, const G4String& pathname, const G4String& filename);

//extracts value from a number*dimension (e.g. 10*mm)
G4double getDoubleValue( G4String paramString );

#endif

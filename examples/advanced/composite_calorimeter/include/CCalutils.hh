///////////////////////////////////////////////////////////////////////////////
// File: CCalutils.hh
// Description: General utilities.
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalutils_hh
#define CCalutils_hh 1

#include "g4std/iostream"
#include "g4std/fstream"
#include "globals.hh"


G4String operator+(const G4String&, const int);
G4String operator+(const G4String&, const double);
// "number " + i = "number i"

G4std::ifstream& readName(G4std::ifstream&, G4String&);
// It reads a name into G4String between quotes and skips lines 
// beginning with #. and if found *ENDDO returns.

G4std::ifstream& findDO(G4std::ifstream&, const G4String&);
// It reads until a *DO str is found.

G4std::ostream& tab(G4std::ostream&);
// It add a tab.

G4std::istream& jump(G4std::istream&);
// It ignores character until the end of line.

bool openGeomFile(G4std::ifstream& is, const G4String& pathname, 
		  const G4String& filename);
// It opens the geometry file, either locally (if it exists) or "remotely".


#endif

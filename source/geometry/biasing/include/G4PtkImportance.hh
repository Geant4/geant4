#ifndef G4PtkImportance_hh
#define G4PtkImportance_hh G4PtkImportance_hh

#include "g4std/map"
#include "globals.hh"
#include "G4PTouchableKey.hh"


typedef map<G4PTouchableKey, G4double, G4PTkComp> G4PtkImportance;

ostream& operator<<(ostream &out, const G4PtkImportance &ptki);


#endif

#ifndef G4PtkImportance_hh
#define G4PtkImportance_hh G4PtkImportance_hh

#include "g4std/map"
#include "globals.hh"
#include "G4PTouchableKey.hh"

typedef G4std::map<G4PTouchableKey, G4double, G4PTkComp> G4PtkImportance;

G4std::ostream& operator<<(G4std::ostream &out, const G4PtkImportance &ptki);


#endif

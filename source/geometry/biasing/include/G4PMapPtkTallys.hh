#ifndef G4PMapPtkTallys_hh 
#define G4PMapPtkTallys_hh G4PMapPtkTallys_hh 

#include "g4std/map"
#include "globals.hh"
#include <iostream>
#include "G4PTouchableKey.hh"

class G4Sigma;

typedef G4std::map<const char *, G4Sigma> G4PMapNameTally;

typedef G4std::map<G4PTouchableKey, G4PMapNameTally, G4PTkComp> G4PMapPtkTallys; 

G4std::ostream& operator<<(G4std::ostream &out, const G4PMapNameTally &tally);
G4std::ostream& operator<<(G4std::ostream &out, const G4PMapPtkTallys &ptktally);

#endif

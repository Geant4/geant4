#ifndef G4PMapPtkTallys_hh 
#define G4PMapPtkTallys_hh G4PMapPtkTallys_hh 

#include "g4std/map"
#include "globals.hh"
#include "G4PTouchableKey.hh"


class G4Sigma;

typedef map<const char *, G4Sigma> G4PMapNameTally;

typedef map<G4PTouchableKey, G4PMapNameTally, G4PTkComp> G4PMapPtkTallys; 

ostream& operator<<(ostream &out, const G4PMapNameTally &tally);
ostream& operator<<(ostream &out, const G4PMapPtkTallys &ptktally);

#endif

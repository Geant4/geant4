#include "G4PMapPtkTallys.hh"
#include "G4PStepStream.hh"
#include "G4Sigma.hh"


ostream& operator<<(ostream &out, const G4PMapNameTally &tallys){
  for (G4PMapNameTally::const_iterator it = tallys.begin();
       it != tallys.end(); it++) {
    out << it->first << "\n";
    out << it->second;
  }
  return out;
}
ostream& operator<<(ostream &out, const G4PMapPtkTallys &tallystore){
  for (G4PMapPtkTallys::const_iterator it = tallystore.begin();
       it != tallystore.end(); it++) {
    out << it->first << "\n";
    out << it->second;
  }
  return out;
}

#include "G4PtkImportance.hh"
#include "G4PStepStream.hh"

ostream& operator<<(ostream &out, const G4PtkImportance &ptki){
  for (G4PtkImportance::const_iterator it = ptki.begin();
       it != ptki.end(); it++) {
    out << it->first << ", importance = ";
    out << it->second << "\n";
  }
  return out;
}


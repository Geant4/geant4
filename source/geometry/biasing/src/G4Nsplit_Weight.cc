#include "G4Nsplit_Weight.hh"

G4std::ostream& operator<<(G4std::ostream &out, const G4Nsplit_Weight  &nw){
  out << "nw.fN = " << nw.fN << ", nw.fW = " << nw.fW;
  return out;
}

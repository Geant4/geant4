#include "G4Nsplit_Weight.hh"

ostream& operator<<(ostream &out, const G4Nsplit_Weight  &nw){
  out << "nw.fN = " << nw.fN << ", nw.fW = " << nw.fW;
  return out;
}

#include "G4Pstring.hh"


G4String str(const int &i){
  char s[200];
  sprintf(s,"%d",i);
  return s;
}
G4String str(const double &d){
  char s[200];
  sprintf(s,"%f",d);
  return s;
}


G4String str(const G4ThreeVector &v) {
  G4String out = "(";
  out += str(v.x());
  out += ", "; 
  out += str(v.y()); 
  out += ", ";
  out += str(v.z());
  out += ")";
  return out;
}



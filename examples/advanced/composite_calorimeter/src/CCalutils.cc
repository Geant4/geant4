///////////////////////////////////////////////////////////////////////////////
// File: CCalutils.cc
// Description: General utilities.
///////////////////////////////////////////////////////////////////////////////
#include "CCalutils.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "g4std/strstream"


G4String operator+(const G4String& str, const int i) {
  int l = str.length() + 15; //How long can an integer be?
  char *cname = new char[l];
  cname[0]='\0';
  G4std::ostrstream os(cname, l);
  os << str << i <<'\0';
  G4String back(cname);  
  delete[] cname;
  return back;
}


G4String operator+(const G4String& str, const double i) {
  int l = str.length() + 15; //How long can an double be?
  char *cname = new char[l];
  cname[0]='\0';
  G4std::ostrstream os(cname, l);
  os << str << i <<'\0';
  G4String back(cname);  
  delete[] cname;
  return back;
}


G4std::ifstream& readName(G4std::ifstream& is, G4String& name){
  is >> name;
  if ( name != "*ENDDO" ) {
    while ( name.find("#.") != G4String::npos ) { // It is a comment. Skip line.
      is.ignore(999,'\n');
      is >> name;
    };
    while ( name.rfind('\"') != name.length()-1 ) {
      G4String other;
      is >> other;
      name += " ";
      name += other;
    };
    name = name.strip(G4String::both, '\"');
  }  
  return is;
}


G4std::ifstream& findDO(G4std::ifstream& is, const G4String& str){
  // Loop until *DO str is found
  G4String firstwd, dowhat;
  dowhat = "";
  while (dowhat != str && is) {
    is >> firstwd;
    while (firstwd != "*DO" && is) {
      is.ignore(999,'\n');
      is >> firstwd;
    };
    is >> dowhat;
  };
  is.ignore(999,'\n');
  return is;
}


G4std::ostream& tab(G4std::ostream& os) {
  os << '\t';
  return os;
}


G4std::istream& jump(G4std::istream& is) {
  char first, second;
  is.ignore(999,'\n');
  do {
    is.get(first);
    second = is.peek();
    if (first == '#' && second =='.') {
      is.ignore(999,'\n');
    }
    else if (first == '\n'); //it was already picked
    else {
      is.putback(first);
      return is;
    }
  }
  while (is);
  return is;
}


bool openGeomFile(G4std::ifstream& is, 
		  const G4String& pathname, const G4String& filename) {
  G4String fullname = pathname+"/"+filename;
  is.open( fullname.c_str() );
  if (!is) {
    G4cerr << "ERROR: Could not open file " << filename << G4endl;
    return false;
  }
  return true;
}

///////////////////////////////////////////////////////////////////////////////
// File: CCalutils.cc
// Description: General utilities.
///////////////////////////////////////////////////////////////////////////////
#include "CCalutils.hh"
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

G4std::istream& readName(G4std::istream& is, G4String& name){
  is >> name;
  if (name!="*ENDDO") {
    while (name.index("#.")==0) { //It is a comment. Skip line.
      is.ignore(999,'\n');
      is >> name;
    };
    while (name.last('\"') != name.length()-1) {
      G4String other;
      is >> other;
      name+= " ";
      name+=other;
    };
    name = name.strip(G4String::both, '\"');
  }  
  return is;
}

G4std::istream& findDO(G4std::istream& is, const G4String& str){
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
  //Check first if the file exists loacally
  is.open(filename);
  if (!is) { //if filename does not exist or is not readable...

    //Try the "remote" loacation for the file    
    G4String fullname = pathname+"/"+filename;
    is.open(fullname);
    if (!is) {
      G4cerr << "ERROR: Could not open file " << filename << G4endl;
      return false;
    }
  }
  
  return true;
}


G4double getDoubleValue( G4String paramString )
{
  G4double vl;
  char unts[30];

  const char* t = paramString;
  G4std::istrstream is((char*)t);
  is >> vl >> unts;
  G4String unt = unts;
  // take away the '*'
  return (vl*G4UnitDefinition::GetValueOf( unt.substr(1,unt.size()) ));
}

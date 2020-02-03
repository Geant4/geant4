//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
///////////////////////////////////////////////////////////////////////////////
// File: CCalutils.cc
// Description: General utilities.
///////////////////////////////////////////////////////////////////////////////
#include "CCalutils.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

#include <sstream>


G4String operator+(const G4String& str, const G4int i) {
  std::ostringstream os;
  os << str << i <<'\0';
  G4String back = os.str();  
  return back;
}


G4String operator+(const G4String& str, const G4double i) {
  std::ostringstream os;
  os << str << i <<'\0';
  G4String back = os.str();  
  return back;
}


std::ifstream& readName(std::ifstream& is, G4String& name){
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


std::ifstream& findDO(std::ifstream& is, const G4String& str){
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


std::ostream& tab(std::ostream& os) {
  os << '\t';
  return os;
}


std::istream& jump(std::istream& is) {
  char first = ' ';
  char second = ' ';
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


G4bool openGeomFile(std::ifstream& is, 
                    const G4String& pathname, const G4String& filename) {
  G4String fullname = pathname+"/"+filename;
  is.open( fullname.c_str() );
  if (!is) {
    G4cerr << "ERROR: Could not open file " << filename << G4endl;
    return false;
  }
  return true;
}

//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4StringConversion.cc,v 1.1 2002-10-16 14:30:00 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4Pstring.cc
//
// ----------------------------------------------------------------------

#include "G4StringConversion.hh"
#include "g4std/strstream"
#include "G4Nsplit_Weight.hh"


G4String str(G4int i){
  //  G4std::ostrstream os;
  //  os << i << "\0";
  //  G4String s(os.str());
  char s[200];
  G4std::sprintf(s,"%d",i);
  return s;
}

G4String str(G4double d){
  //  G4std::ostrstream os;
  //  os << d << "\0";
  //  G4String s(os.str());
  char s[200];
  G4std::sprintf(s,"%f",d);
  return s;
}

G4String str(const G4ThreeVector &v){
  //  G4std::ostrstream os;
  //  os << "(" 
  //     << G4std::str(v.x()) 
  //     << ", "  
  //     << G4std::str(v.y())
  //     << ", "
  //     << G4std::str(v.z())
  //     << ")\0";
  //  G4String s(os.str());
  G4String out = G4String("(");
  out += G4std::str(v.x());
  out += ", "; 
  out += G4std::str(v.y()); 
  out += ", ";
  out += G4std::str(v.z());
  out += ")";
  return out;
}

G4std::ostream& operator<<(G4std::ostream &out, const G4Nsplit_Weight &nw)
{
  out << "nw.fN = " << nw.fN << ", nw.fW = " << nw.fW;
  return out;
}

G4String str(const G4Nsplit_Weight &nw){
  G4String out = G4String("nw.fN = ");
  out += G4std::str(nw.fN);
  out += ", nw.fW = "; 
  out += G4std::str(nw.fW); 
  return out;
 
}



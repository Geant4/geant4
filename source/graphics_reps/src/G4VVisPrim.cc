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
// $Id: G4VVisPrim.cc,v 1.7 2001-07-11 10:01:09 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  August 1995

#include "G4VVisPrim.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

G4VVisPrim::G4VVisPrim () {}

G4VVisPrim::G4VVisPrim (const G4VVisPrim& prim):
G4Visible (prim)
{}

G4VVisPrim::G4VVisPrim (const G4VisAttributes* pVA):
G4Visible (pVA)
{}

G4VVisPrim::~G4VVisPrim () {}

G4Visible& G4VVisPrim::operator = (const G4Visible& right) {
  return G4Visible::operator = (right);
}

G4VVisPrim& G4VVisPrim::operator = (const G4VVisPrim& right) {
  if (&right == this) return *this;
  G4Visible::operator = (right);
  return *this;
}

G4bool G4VVisPrim::operator == (const G4Visible& right) const{
  return G4Visible::operator == (right);
}

G4bool G4VVisPrim::operator == (const G4VVisPrim& right) const{
  return G4Visible::operator == (right);
}

G4std::ostream& operator << (G4std::ostream& os, const G4VVisPrim& prim) {
  return os << (G4Visible) prim;
}

// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VVisPrim.cc,v 2.1 1998/07/13 16:56:21 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// John Allison  August 1995

#include "G4VVisPrim.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

ostream& operator << (ostream& os, const G4VVisPrim& prim) {
  return os << (G4Visible) prim;
}

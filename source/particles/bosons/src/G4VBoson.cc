// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VBoson.cc,v 1.1 1999-01-07 16:09:49 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
// --------------------------------------------------------------

#include "G4VBoson.hh"

const G4VBoson & G4VBoson::operator=(const G4VBoson &right)
{
  if (this != &right) {
  } return right;
}

G4int G4VBoson::operator==(const G4VBoson &right) const
{
  return (this->GetParticleName() == right.GetParticleName());
}

G4int G4VBoson::operator!=(const G4VBoson &right) const
{
  return (this->GetParticleName() != right.GetParticleName());
}

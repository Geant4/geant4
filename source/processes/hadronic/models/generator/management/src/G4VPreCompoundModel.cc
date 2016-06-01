// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VPreCompoundModel.cc,v 1.3 1998/09/28 09:25:23 hpw Exp $
// GEANT4 tag $Name: geant4-00 $
//

#include "G4VPreCompoundModel.hh"

G4VPreCompoundModel::G4VPreCompoundModel(G4ExcitationHandler *const value):
  theExcitationHandler(value)
{
}


const G4VPreCompoundModel & 
G4VPreCompoundModel::operator=(const G4VPreCompoundModel &right)
{
  G4Exception("G4VPreCompoundModel::operator= meant to not be accessable");
  return *this;
}

G4bool G4VPreCompoundModel::operator==(const G4VPreCompoundModel &right) const
{
  return false;
}

G4bool G4VPreCompoundModel::operator!=(const G4VPreCompoundModel &right) const
{
  return true;
}

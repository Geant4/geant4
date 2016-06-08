// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VPreCompoundModel.cc,v 1.2 1999/12/15 14:52:38 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
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

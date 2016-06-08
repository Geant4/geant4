// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4StringModel.cc,v 1.1.10.1 1999/12/07 20:51:54 gunter Exp $
// GEANT4 tag $Name: geant4-01-01 $
//
// G4StringModel
#include "G4StringModel.hh"

G4StringModel::G4StringModel()
{
	the3DNucleus=NULL;
	theStringFragmentationModel=NULL;
	theGenerator=NULL;
}

G4StringModel::G4StringModel(const G4StringModel &right)
{
}


G4StringModel::~G4StringModel()
{
}


const G4StringModel & G4StringModel::operator=(const G4StringModel &right)
{
  G4Exception("G4StringModel::operator= meant to not be accessable");
  return *this;
}


int G4StringModel::operator==(const G4StringModel &right) const
{
  return 0;
}

int G4StringModel::operator!=(const G4StringModel &right) const
{
  return 1;
}

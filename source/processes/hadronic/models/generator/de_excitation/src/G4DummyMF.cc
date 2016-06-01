// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DummyMF.cc,v 1.1 1998/08/22 08:53:45 hpw Exp $
// GEANT4 tag $Name: geant4-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (May 1998)

#include "G4DummyMF.hh"

G4DummyMF::G4DummyMF()
{
}

G4DummyMF::G4DummyMF(const G4DummyMF &right)
{
}

G4DummyMF::~G4DummyMF()
{
}

const G4DummyMF & G4DummyMF::operator=(const G4DummyMF &right)
{
  G4Exception("G4DummyMF::operator= meant to not be accessable");
  return *this;
}

int G4DummyMF::operator==(const G4DummyMF &right) const
{
  return 0;
}

int G4DummyMF::operator!=(const G4DummyMF &right) const
{
  return 1;
}

G4FragmentVector * G4DummyMF::BreakItUp(const G4Fragment &theNucleus)
{
//  G4cout << "G4DummyMF::BreakItUp called"<<endl;
  G4FragmentVector * theResult =
           new G4FragmentVector;
  // all calculations here
  return theResult;
}

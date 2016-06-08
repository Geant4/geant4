// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EventGenerator.cc,v 1.1.10.1 1999/12/07 20:51:53 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
// G4EventGenerator
#include "G4EventGenerator.hh"

G4EventGenerator::G4EventGenerator()
{
   SetMinEnergy (0 *GeV);
   SetMaxEnergy (0 *GeV);
}

G4EventGenerator::G4EventGenerator(const G4EventGenerator &right)
{
}


G4EventGenerator::~G4EventGenerator()
{
}


const G4EventGenerator & G4EventGenerator::operator=(const G4EventGenerator &right)
{
  G4Exception("G4EventGenerator::operator= meant to not be accessable");
  return *this;
}


int G4EventGenerator::operator==(const G4EventGenerator &right) const
{
  return 0;
}

int G4EventGenerator::operator!=(const G4EventGenerator &right) const
{
  return 1;
}

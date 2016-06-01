// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VHighEnergyGenerator.cc,v 1.1 1998/08/22 08:55:56 hpw Exp $
// GEANT4 tag $Name: geant4-00 $
//
// G4VHighEnergyGenerator
#include "G4VHighEnergyGenerator.hh"

G4VHighEnergyGenerator::G4VHighEnergyGenerator()
{
}

G4VHighEnergyGenerator::G4VHighEnergyGenerator(const G4VHighEnergyGenerator &right)
{
}


G4VHighEnergyGenerator::~G4VHighEnergyGenerator()
{
}


const G4VHighEnergyGenerator & G4VHighEnergyGenerator::operator=(const G4VHighEnergyGenerator &right)
{
  G4Exception("G4VHighEnergyGenerator::operator= meant to not be accessable");
  return *this;
}


int G4VHighEnergyGenerator::operator==(const G4VHighEnergyGenerator &right) const
{
  return 0;
}

int G4VHighEnergyGenerator::operator!=(const G4VHighEnergyGenerator &right) const
{
  return 1;
}


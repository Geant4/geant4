// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VKineticNucleon.cc,v 1.1 1999/04/15 12:10:12 hpw Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
#include "G4VKineticNucleon.hh"

G4VKineticNucleon::G4VKineticNucleon()
{
}

G4VKineticNucleon::G4VKineticNucleon(const G4VKineticNucleon &right)
{
}


G4VKineticNucleon::~G4VKineticNucleon()
{
}


//const G4VKineticNucleon & G4VKineticNucleon::operator=(const G4VKineticNucleon &right)
//{}


int G4VKineticNucleon::operator==(const G4VKineticNucleon &right) const
{
	return this == &right;
}

int G4VKineticNucleon::operator!=(const G4VKineticNucleon &right) const
{
	return this != &right;

}





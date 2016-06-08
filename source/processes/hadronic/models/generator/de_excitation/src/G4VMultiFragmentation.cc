// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#include "G4VMultiFragmentation.hh"


G4VMultiFragmentation::G4VMultiFragmentation()
{
}

G4VMultiFragmentation::G4VMultiFragmentation(const G4VMultiFragmentation &right)
{
	G4Exception("G4VMultiFragmentation::copy_constructor meant to not be accessable");
}


G4VMultiFragmentation::~G4VMultiFragmentation()
{
}


const G4VMultiFragmentation & G4VMultiFragmentation::operator=(const G4VMultiFragmentation &right)
{
	G4Exception("G4VMultiFragmentation::operator= meant to not be accessable");
  	return *this;
}


G4bool G4VMultiFragmentation::operator==(const G4VMultiFragmentation &right) const
{
	G4Exception("G4VMultiFragmentation::operator== meant to not be accessable");
	return false;
}

G4bool G4VMultiFragmentation::operator!=(const G4VMultiFragmentation &right) const
{
	G4Exception("G4VMultiFragmentation::operator=! meant to not be accessable");
	return true;
}




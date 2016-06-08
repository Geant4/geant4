// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov. 1999)
//

#include "G4AlphaEvaporationChannel.hh"


const G4AlphaEvaporationChannel & G4AlphaEvaporationChannel::operator=(const G4AlphaEvaporationChannel & right)
{
	G4Exception("G4AlphaEvaporationChannel::operator= meant to not be accessable");
	return *this;
}

G4AlphaEvaporationChannel::G4AlphaEvaporationChannel(const G4AlphaEvaporationChannel & right)
{
	G4Exception("G4AlphaEvaporationChannel::CopyConstructor meant to not be accessable");
}

G4bool G4AlphaEvaporationChannel::operator==(const G4AlphaEvaporationChannel & right) const 
{
	return (this == (G4AlphaEvaporationChannel *) &right);
	//  return false;
}

G4bool G4AlphaEvaporationChannel::operator!=(const G4AlphaEvaporationChannel & right) const 
{
	return (this != (G4AlphaEvaporationChannel *) &right);
	//  return true;
}


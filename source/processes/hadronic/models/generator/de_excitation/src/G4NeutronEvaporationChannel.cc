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

#include "G4NeutronEvaporationChannel.hh"


const G4NeutronEvaporationChannel & G4NeutronEvaporationChannel::operator=(const G4NeutronEvaporationChannel & right)
{
	G4Exception("G4NeutronEvaporationChannel::operator= meant to not be accessable");
	return *this;
}

G4NeutronEvaporationChannel::G4NeutronEvaporationChannel(const G4NeutronEvaporationChannel & right)
{
	G4Exception("G4NeutronEvaporationChannel::CopyConstructor meant to not be accessable");
}

G4bool G4NeutronEvaporationChannel::operator==(const G4NeutronEvaporationChannel & right) const 
{
	return (this == (G4NeutronEvaporationChannel *) &right);
	//  return false;
}

G4bool G4NeutronEvaporationChannel::operator!=(const G4NeutronEvaporationChannel & right) const 
{
	return (this != (G4NeutronEvaporationChannel *) &right);
	//  return true;
}


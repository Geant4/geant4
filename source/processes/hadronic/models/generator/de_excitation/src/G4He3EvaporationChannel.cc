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

#include "G4He3EvaporationChannel.hh"


const G4He3EvaporationChannel & G4He3EvaporationChannel::operator=(const G4He3EvaporationChannel & right)
{
	G4Exception("G4He3EvaporationChannel::operator= meant to not be accessable");
	return *this;
}

G4He3EvaporationChannel::G4He3EvaporationChannel(const G4He3EvaporationChannel & right)
{
	G4Exception("G4He3EvaporationChannel::CopyConstructor meant to not be accessable");
}

G4bool G4He3EvaporationChannel::operator==(const G4He3EvaporationChannel & right) const 
{
	return (this == (G4He3EvaporationChannel *) &right);
	//  return false;
}

G4bool G4He3EvaporationChannel::operator!=(const G4He3EvaporationChannel & right) const 
{
	return (this != (G4He3EvaporationChannel *) &right);
	//  return true;
}


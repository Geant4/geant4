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

#include "G4ProtonEvaporationChannel.hh"


const G4ProtonEvaporationChannel & G4ProtonEvaporationChannel::operator=(const G4ProtonEvaporationChannel & right)
{
	G4Exception("G4ProtonEvaporationChannel::operator= meant to not be accessable");
	return *this;
}

G4ProtonEvaporationChannel::G4ProtonEvaporationChannel(const G4ProtonEvaporationChannel & right)
{
	G4Exception("G4ProtonEvaporationChannel::CopyConstructor meant to not be accessable");
}

G4bool G4ProtonEvaporationChannel::operator==(const G4ProtonEvaporationChannel & right) const 
{
	return (this == (G4ProtonEvaporationChannel *) &right);
	//  return false;
}

G4bool G4ProtonEvaporationChannel::operator!=(const G4ProtonEvaporationChannel & right) const 
{
	return (this != (G4ProtonEvaporationChannel *) &right);
	//  return true;
}



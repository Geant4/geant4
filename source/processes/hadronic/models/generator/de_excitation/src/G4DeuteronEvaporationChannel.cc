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

#include "G4DeuteronEvaporationChannel.hh"


const G4DeuteronEvaporationChannel & G4DeuteronEvaporationChannel::operator=(const G4DeuteronEvaporationChannel & right)
{
	G4Exception("G4DeuteronEvaporationChannel::operator= meant to not be accessable");
	return *this;
}

G4DeuteronEvaporationChannel::G4DeuteronEvaporationChannel(const G4DeuteronEvaporationChannel & right)
{
	G4Exception("G4DeuteronEvaporationChannel::CopyConstructor meant to not be accessable");
}

G4bool G4DeuteronEvaporationChannel::operator==(const G4DeuteronEvaporationChannel & right) const 
{
	return (this == (G4DeuteronEvaporationChannel *) &right);
	//  return false;
}

G4bool G4DeuteronEvaporationChannel::operator!=(const G4DeuteronEvaporationChannel & right) const 
{
	return (this != (G4DeuteronEvaporationChannel *) &right);
	//  return true;
}

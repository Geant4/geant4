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

#include "G4TritonEvaporationChannel.hh"


const G4TritonEvaporationChannel & G4TritonEvaporationChannel::operator=(const G4TritonEvaporationChannel & right)
{
	G4Exception("G4TritonEvaporationChannel::operator= meant to not be accessable");
	return *this;
}

G4TritonEvaporationChannel::G4TritonEvaporationChannel(const G4TritonEvaporationChannel & right)
{
	G4Exception("G4TritonEvaporationChannel::CopyConstructor meant to not be accessable");
}

G4bool G4TritonEvaporationChannel::operator==(const G4TritonEvaporationChannel & right) const 
{
	return (this == (G4TritonEvaporationChannel *) &right);
	//  return false;
}

G4bool G4TritonEvaporationChannel::operator!=(const G4TritonEvaporationChannel & right) const 
{
	return (this != (G4TritonEvaporationChannel *) &right);
	//  return true;
}

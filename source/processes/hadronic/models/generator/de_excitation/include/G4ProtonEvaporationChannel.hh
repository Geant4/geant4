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


#ifndef G4ProtonEvaporationChannel_h
#define G4ProtonEvaporationChannel_h 1

#include "G4EvaporationChannel.hh"
#include "G4ProtonCoulombBarrier.hh"
#include "G4ProtonEvaporationProbability.hh"

class G4ProtonEvaporationChannel : public G4EvaporationChannel
{
public:
	// only available constructor
	G4ProtonEvaporationChannel() : G4EvaporationChannel(1,1,
		&theEvaporationProbability,&theCoulombBarrier) {};

	// destructor
	~G4ProtonEvaporationChannel() {};

private:
	const G4ProtonEvaporationChannel & operator=(const G4ProtonEvaporationChannel & right);  

	G4ProtonEvaporationChannel(const G4ProtonEvaporationChannel & right);

public:
	G4bool operator==(const G4ProtonEvaporationChannel & right) const;
	G4bool operator!=(const G4ProtonEvaporationChannel & right) const;


private:

	G4ProtonEvaporationProbability  theEvaporationProbability;
	
	G4ProtonCoulombBarrier  theCoulombBarrier;
};
#endif

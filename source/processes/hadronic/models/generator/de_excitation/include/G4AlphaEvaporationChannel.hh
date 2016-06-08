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


#ifndef G4AlphaEvaporationChannel_h
#define G4AlphaEvaporationChannel_h 1

#include "G4EvaporationChannel.hh"
#include "G4AlphaCoulombBarrier.hh"
#include "G4AlphaEvaporationProbability.hh"

class G4AlphaEvaporationChannel : public G4EvaporationChannel
{
public:
	// only available constructor
	G4AlphaEvaporationChannel() : G4EvaporationChannel(4,2,
		&theEvaporationProbability,&theCoulombBarrier) {};

	// destructor
	~G4AlphaEvaporationChannel() {};

private:
	const G4AlphaEvaporationChannel & operator=(const G4AlphaEvaporationChannel & right);  

	G4AlphaEvaporationChannel(const G4AlphaEvaporationChannel & right);


public:
	G4bool operator==(const G4AlphaEvaporationChannel & right) const;
	G4bool operator!=(const G4AlphaEvaporationChannel & right) const;

private:

	G4AlphaCoulombBarrier theCoulombBarrier;
	
	G4AlphaEvaporationProbability theEvaporationProbability;

};
#endif

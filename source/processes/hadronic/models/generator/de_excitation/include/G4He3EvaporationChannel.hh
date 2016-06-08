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


#ifndef G4He3EvaporationChannel_h
#define G4He3EvaporationChannel_h 1

#include "G4EvaporationChannel.hh"
#include "G4He3CoulombBarrier.hh"
#include "G4He3EvaporationProbability.hh"

class G4He3EvaporationChannel : public G4EvaporationChannel
{
public:
	// only available constructor
	G4He3EvaporationChannel() : G4EvaporationChannel(3,2,
		&theEvaporationProbability,&theCoulombBarrier) {};

	// destructor
	~G4He3EvaporationChannel() {};

private:
	const G4He3EvaporationChannel & operator=(const G4He3EvaporationChannel & right);  

	G4He3EvaporationChannel(const G4He3EvaporationChannel & right);

public:
	G4bool operator==(const G4He3EvaporationChannel & right) const;
	G4bool operator!=(const G4He3EvaporationChannel & right) const;

private:

	G4He3CoulombBarrier theCoulombBarrier;
	
	G4He3EvaporationProbability theEvaporationProbability;

};
#endif

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


#ifndef G4NeutronEvaporationChannel_h
#define G4NeutronEvaporationChannel_h 1

#include "G4EvaporationChannel.hh"
#include "G4NeutronCoulombBarrier.hh"
#include "G4NeutronEvaporationProbability.hh"

class G4NeutronEvaporationChannel : public G4EvaporationChannel
{
public:
	// only available constructor
	G4NeutronEvaporationChannel() : G4EvaporationChannel(1,0,
		&theEvaporationProbability,&theCoulombBarrier) {};

	// destructor
	~G4NeutronEvaporationChannel() {};

private:
	const G4NeutronEvaporationChannel & operator=(const G4NeutronEvaporationChannel & right);  

	G4NeutronEvaporationChannel(const G4NeutronEvaporationChannel & right);

public:
	G4bool operator==(const G4NeutronEvaporationChannel & right) const;
	G4bool operator!=(const G4NeutronEvaporationChannel & right) const;

private:

	G4NeutronCoulombBarrier theCoulombBarrier;
	
	G4NeutronEvaporationProbability theEvaporationProbability;

};
#endif

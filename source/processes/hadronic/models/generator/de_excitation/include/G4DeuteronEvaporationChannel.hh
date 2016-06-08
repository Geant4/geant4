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


#ifndef G4DeuteronEvaporationChannel_h
#define G4DeuteronEvaporationChannel_h 1

#include "G4EvaporationChannel.hh"
#include "G4DeuteronCoulombBarrier.hh"
#include "G4DeuteronEvaporationProbability.hh"

class G4DeuteronEvaporationChannel : public G4EvaporationChannel
{
public:
	// only available constructor
	G4DeuteronEvaporationChannel() : G4EvaporationChannel(2,1,
		&theEvaporationProbability,&theCoulombBarrier) {};

	// destructor
	~G4DeuteronEvaporationChannel() {};

private:
	const G4DeuteronEvaporationChannel & operator=(const G4DeuteronEvaporationChannel & right);  

	G4DeuteronEvaporationChannel(const G4DeuteronEvaporationChannel & right);

public:
	G4bool operator==(const G4DeuteronEvaporationChannel & right) const;
	G4bool operator!=(const G4DeuteronEvaporationChannel & right) const;

private:

	G4DeuteronCoulombBarrier theCoulombBarrier;
	
	G4DeuteronEvaporationProbability theEvaporationProbability;

};
#endif

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


#ifndef G4TritonEvaporationChannel_h
#define G4TritonEvaporationChannel_h 1

#include "G4EvaporationChannel.hh"
#include "G4TritonCoulombBarrier.hh"
#include "G4TritonEvaporationProbability.hh"

class G4TritonEvaporationChannel : public G4EvaporationChannel
{
public:
	// only available constructor
	G4TritonEvaporationChannel() : G4EvaporationChannel(3,1,
		&theEvaporationProbability,&theCoulombBarrier) {};

	// destructor
	~G4TritonEvaporationChannel() {};

private:
	const G4TritonEvaporationChannel & operator=(const G4TritonEvaporationChannel & right);  

	G4TritonEvaporationChannel(const G4TritonEvaporationChannel & right);

public:
	G4bool operator==(const G4TritonEvaporationChannel & right) const;
	G4bool operator!=(const G4TritonEvaporationChannel & right) const;

private:

	G4TritonCoulombBarrier theCoulombBarrier;
	
	G4TritonEvaporationProbability theEvaporationProbability;

};
#endif

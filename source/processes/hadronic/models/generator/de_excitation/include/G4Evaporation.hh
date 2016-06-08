// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations

#ifndef G4Evaporation_h
#define G4Evaporation_h 1

#include "globals.hh"
#include "g4rw/tvvector.h"
#include "g4rw/tpordvec.h"
#include "g4rw/tvordvec.h"

#include "G4ios.hh"
#include "G4VEvaporation.hh"
#include "G4VEvaporationChannel.hh"
#include "G4Fragment.hh"
#include "G4NucleiProperties.hh"
#include "Randomize.hh"


class G4Evaporation : public G4VEvaporation
{
public:
	G4Evaporation();
	G4Evaporation(G4RWTPtrOrderedVector<G4VEvaporationChannel> * aChannelsVector) :
		theChannels(aChannelsVector),
		myOwnChannelsVector(false)
	 {};
	 
	~G4Evaporation();

private:
	G4Evaporation(const G4Evaporation &right);

	const G4Evaporation & operator=(const G4Evaporation &right);
	G4bool operator==(const G4Evaporation &right) const;
	G4bool operator!=(const G4Evaporation &right) const;

public:
	G4FragmentVector * BreakItUp(const G4Fragment &theNucleus);
		
private:

	G4bool myOwnChannelsVector;
	G4RWTPtrOrderedVector<G4VEvaporationChannel> * theChannels;
};

#endif






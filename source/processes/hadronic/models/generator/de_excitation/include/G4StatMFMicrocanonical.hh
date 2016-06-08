#ifndef G4StatMFMicroCanonical_h
#define G4StatMFMicroCanonical_h 1

//#include <rw/tvvector.h>
#include "g4rw/tvordvec.h"

#include "G4VStatMFEnsemble.hh"
#include "G4StatMFMicroPartition.hh"
#include "G4StatMFMicroManager.hh"
#include "G4StatMFParameters.hh"
#include "G4StatMFChannel.hh"

#include "G4Fragment.hh"
#include "Randomize.hh"


class G4StatMFMicroCanonical : public G4VStatMFEnsemble {

public:

	// G4StatMFMicroCanonical class must be initialized with a G4Fragment.
	G4StatMFMicroCanonical(const G4Fragment & theFragment);

	// destructor
	~G4StatMFMicroCanonical();

private:
	// default constructor
	G4StatMFMicroCanonical() {};


	// copy constructor
	G4StatMFMicroCanonical(const G4StatMFMicroCanonical &right);


	// operators
	G4StatMFMicroCanonical & operator=(const G4StatMFMicroCanonical & right);
	G4bool operator==(const G4StatMFMicroCanonical & right) const;
	G4bool operator!=(const G4StatMFMicroCanonical & right) const;


public:

	// Choice of fragment atomic numbers and charges.
	G4StatMFChannel * ChooseAandZ(const G4Fragment & theFragment);
	
	enum {MaxAllowedMultiplicity = 4};

private:

	// Initailization method
	void Initialize(const G4Fragment & theFragment);

	// Calculate Entropy of Compound Nucleus
	G4double CalcEntropyOfCompoundNucleus(const G4Fragment & theFragment, G4double & TConf);


	G4double CalcFreeInternalEnergy(const G4Fragment & theFragment, const G4double T);


	G4double CalcInvLevelDensity(const G4int anA);
	
	
// Data members
private:
	
	// This is a vector of partitions managers for partitions of different 
	// multiplicities:
	
	G4RWTPtrOrderedVector<G4StatMFMicroManager> _ThePartitionManagerVector;
	


	// Statistical weight of compound nucleus
	G4double _WCompoundNucleus;

};

#endif

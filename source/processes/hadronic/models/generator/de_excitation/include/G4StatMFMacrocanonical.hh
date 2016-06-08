#ifndef G4StatMFMacroCanonical_h
#define G4StatMFMacroCanonical_h 1

#include "g4rw/tpordvec.h"

#include "G4Fragment.hh"
#include "G4StatMFFragment.hh"
#include "G4VStatMFEnsemble.hh"
#include "G4VStatMFMacroCluster.hh"
#include "G4StatMFMacroNucleon.hh"
#include "G4StatMFMacroBiNucleon.hh"
#include "G4StatMFMacroTriNucleon.hh"
#include "G4StatMFMacroTetraNucleon.hh"
#include "G4StatMFMacroMultiNucleon.hh"
#include "G4StatMFParameters.hh"
#include "G4StatMFChannel.hh"
#include "G4StatMFMacroTemperature.hh"
#include "Randomize.hh"


class G4StatMFMacroCanonical : public G4VStatMFEnsemble {

public:

	// G4StatMFMacroCanonical class must be initialized with a G4Fragment.
	G4StatMFMacroCanonical(const G4Fragment & theFragment);

	// destructor
	~G4StatMFMacroCanonical();

private:
	// default constructor
	G4StatMFMacroCanonical() {};


	// copy constructor
	G4StatMFMacroCanonical(const G4StatMFMacroCanonical &right) {};


	// operators
	G4StatMFMacroCanonical & operator=(const G4StatMFMacroCanonical & right);
	G4bool operator==(const G4StatMFMacroCanonical & right) const;
	G4bool operator!=(const G4StatMFMacroCanonical & right) const;


public:

	// Choice of fragment atomic numbers and charges.
	G4StatMFChannel * ChooseAandZ(const G4Fragment &theFragment);

private:

	// Initailization method
	void Initialize(const G4Fragment & theFragment);

	//
	void CalculateTemperature(const G4Fragment & theFragment);

	// Determines fragments multiplicities and compute total fragment multiplicity
	G4double ChooseA(const G4double A, G4RWTValVector<G4double> & ANumbers);
	
	// Samples charges of fragments
	G4StatMFChannel * ChooseZ(const G4int & Z, 
									G4RWTValOrderedVector<G4double> & FragmentsA);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


	// Chemical Potential \mu
	G4double _ChemPotentialMu;

	// Chemical Potential \nu
	G4double _ChemPotentialNu;


	// Parameter Kappa
	G4double _Kappa;

	// Clusters
	G4RWTPtrOrderedVector<G4VStatMFMacroCluster> _theClusters;

};

#endif

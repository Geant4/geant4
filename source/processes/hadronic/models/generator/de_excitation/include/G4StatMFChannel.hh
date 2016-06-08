#ifndef G4StatMFChannel_h
#define G4StatMFChannel_h 1

#include "G4StatMFParameters.hh"
#include "G4StatMFFragment.hh"


class G4StatMFChannel {

public:
	// Default Constructor
	G4StatMFChannel() : 
	_NumOfNeutralFragments(0), 
	_NumOfChargedFragments(0)
	{}


	// Destructor
	~G4StatMFChannel() { _theFragments.clearAndDestroy(); }


private:

	// Copy constructor
	G4StatMFChannel(const G4StatMFChannel & right);

	// operators
	G4StatMFChannel & operator=(const G4StatMFChannel & right);

public:
	G4bool operator==(const G4StatMFChannel & right) const;
	G4bool operator!=(const G4StatMFChannel & right) const;
	
public:

	void CreateFragment(const G4double A, const G4double Z);
// 	{
// 		_theFragments.insert(new G4StatMFFragment(A,Z));
// 	}
	
	G4int GetMultiplicity(void) { return _theFragments.entries();}
	
	// Return false if there is some unphysical fragment
	G4bool CheckFragments(void);


	G4double GetFragmentsCoulombEnergy(void);


	G4double GetFragmentsEnergy(const G4double T) const;
	
	
	G4FragmentVector * GetFragments(const G4double anA, const G4double anZ, const G4double T);
	

private:


	// This method calculates asymptotic fragments momenta.
	void CoulombImpulse(const G4double anA, const G4double anZ, const G4double T);
	
	void PlaceFragments(const G4double anA);

	void SolveEqOfMotion(const G4double anA, const G4double anZ, const G4double T);


	// Calculates fragments momentum components at the breakup instant.
	// Fragment kinetic energies will be calculated according to the
	// Boltzamann distribution at given temperature.
	void FragmentsMomenta(const G4int NF, const G4int idx, const G4double T);	


	// Samples a isotropic random vectorwith a magnitud given by Magnitude.
	// By default Magnitude = 1
	G4ThreeVector IsotropicVector(const G4double Magnitude = 1.0);


	// Rotates a 3-vector P to close momentum triangle Pa + V + P = 0
	G4ThreeVector RotateMomentum(G4ThreeVector Pa, G4ThreeVector V, 
										  G4ThreeVector P);

private:

	G4RWTPtrOrderedVector<G4StatMFFragment> _theFragments;

	G4int _NumOfNeutralFragments;
	
	G4int _NumOfChargedFragments;

};

#endif

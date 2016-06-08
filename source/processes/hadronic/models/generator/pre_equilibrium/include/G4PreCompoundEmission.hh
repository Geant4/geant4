
#ifndef G4PreCompoundEmission_h
#define G4PreCompoundEmission_h 1

#include "G4VPreCompoundFragment.hh"
#include "G4PreCompoundFragmentVector.hh"
#include "G4ReactionProduct.hh"
#include "G4Fragment.hh"
#include "Randomize.hh"

class G4PreCompoundEmission
{
public:
	G4PreCompoundEmission() {};
	~G4PreCompoundEmission() {};

private:
  G4PreCompoundEmission(const G4PreCompoundEmission &right);
  const G4PreCompoundEmission& operator=(const G4PreCompoundEmission &right);
  G4bool operator==(const G4PreCompoundEmission &right) const;
  G4bool operator!=(const G4PreCompoundEmission &right) const;

public:

	void Initialize(const G4Fragment & aFragment) 
	{
		theFragmentsVector.Initialize(aFragment);
		return;
	}
	
	G4double GetTotalProbability(const G4Fragment & aFragment) 
	{
		return theFragmentsVector.CalculateProbabilities(aFragment);
	}
	
	G4ReactionProduct * PerformEmission(G4Fragment & aFragment);

	
private:

	G4ThreeVector IsotropicRandom3Vector(G4double Magnitude = 1.0) const;

	G4ParticleMomentum RotateMomentum(G4ParticleMomentum Pa, G4ParticleMomentum V,
												 G4ParticleMomentum P) const;


	// A vector with the allowed emission fragments 
	G4PreCompoundFragmentVector theFragmentsVector;

};
#endif

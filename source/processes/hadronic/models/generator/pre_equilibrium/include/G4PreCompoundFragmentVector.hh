#ifndef G4PreCompoundFragmentVector_h
#define G4PreCompoundFragmentVector_h 1


#include "G4VPreCompoundFragment.hh"
#include "g4rw/tpordvec.h"

class G4PreCompoundFragmentVector 
{
public:
	G4PreCompoundFragmentVector();
	~G4PreCompoundFragmentVector();
	
private:
	G4PreCompoundFragmentVector(const G4PreCompoundFragmentVector &right);
	const G4PreCompoundFragmentVector& operator=(const G4PreCompoundFragmentVector &right);
	G4bool operator==(const G4PreCompoundFragmentVector &right) const;
	G4bool operator!=(const G4PreCompoundFragmentVector &right) const;	

public:

	void Initialize(const G4Fragment & aFragment)
	{
		TotalEmissionProbability = 0.0;
		for (G4int i=0; i < theChannels.entries(); i++)	theChannels(i)->Init(aFragment);
		return;
	}
	
	G4double CalculateProbabilities(const G4Fragment & aFragment);
	
	G4VPreCompoundFragment * ChooseFragment(void);
		
private:

	G4RWTPtrOrderedVector<G4VPreCompoundFragment> theChannels;

	G4double TotalEmissionProbability;

};
#endif

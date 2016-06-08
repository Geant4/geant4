#ifndef G4VStatMFEnsemble_h
#define G4VStatMFEnsemble_h 1

#include "G4StatMFParameters.hh"
#include "G4StatMFChannel.hh"

class G4VStatMFEnsemble {

public:
	// Default Constructor
	G4VStatMFEnsemble() :
	__FreeInternalE0(0.0),
	__MeanTemperature(0.0),
	__MeanEntropy(0.0),
	__MeanMultiplicity(0.0)
	{};


	// Destructor
	virtual ~G4VStatMFEnsemble() {};


private:

	// Copy constructor
	G4VStatMFEnsemble(const G4VStatMFEnsemble & right);

	// operators
	G4VStatMFEnsemble & operator=(const G4VStatMFEnsemble & right);
	G4bool operator==(const G4VStatMFEnsemble & right) const;
	G4bool operator!=(const G4VStatMFEnsemble & right) const;
	
public:

	virtual G4StatMFChannel * ChooseAandZ(const G4Fragment & aFragment) = 0;
		
	G4double GetMeanMultiplicity(void) const {return __MeanMultiplicity;}
	
	G4double GetMeanTemperature(void) const {return __MeanTemperature;}

protected:

	// Free internal energy at temperature T = 0
	G4double __FreeInternalE0;
	
	// Mean temperature 
	G4double __MeanTemperature;
	
	// Mean Entropy 
	G4double __MeanEntropy;
	
	// Mean Multiplicity
	G4double __MeanMultiplicity;
};

#endif

// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998) 
//



#ifndef G4EvaporationProbability_h
#define G4EvaporationProbability_h 1


#include "G4VEmissionProbability.hh"
#include "G4VLevelDensityParameter.hh"
#include "G4EvaporationLevelDensityParameter.hh"

#include "g4rw/tvvector.h"

class G4EvaporationProbability : public G4VEmissionProbability
{
public:
	// Only available constructor
	G4EvaporationProbability(const G4int anA, const G4int aZ, const G4double aGamma) : 
		theA(anA),
		theZ(aZ),
		Gamma(aGamma) 
		{
			theEvapLDPptr = new G4EvaporationLevelDensityParameter;
		}

	~G4EvaporationProbability() 
	{
		if (theEvapLDPptr != 0) delete theEvapLDPptr;
	}


	
	G4double GetZ(void) const { return theZ; }
	
	G4double GetA(void) const { return theA;} 

protected:

	void SetExcitationEnergiesPtr(G4RWTValVector<G4double> * anExcitationEnergiesPtr) 
		{ExcitationEnergies = anExcitationEnergiesPtr;}

	void SetExcitationSpinsPtr(G4RWTValVector<G4int> * anExcitationSpinsPtr)
		{ExcitationSpins = anExcitationSpinsPtr;}

  
	// Default constructor
	G4EvaporationProbability() {}
private:
	// Copy constructor
	G4EvaporationProbability(const G4EvaporationProbability &right);

	const G4EvaporationProbability & operator=(const G4EvaporationProbability &right);
	G4bool operator==(const G4EvaporationProbability &right) const;
	G4bool operator!=(const G4EvaporationProbability &right) const;
  
public:
	G4double EmissionProbability(const G4Fragment & fragment, const G4double anEnergy);

private:

	G4double CalcProbability(const G4Fragment & fragment, const G4double MaximalKineticEnergy);
	G4double PairingCorrection(const G4int A, const G4int Z) const;
	virtual G4double CCoeficient(const G4double aZ) const {return 0.0;};

	virtual G4double CalcAlphaParam(const G4Fragment & fragment) const {return 1.0;}
	virtual G4double CalcBetaParam(const G4Fragment & fragment) const {return 1.0;}

	// Data Members

	G4VLevelDensityParameter * theEvapLDPptr;
	
	G4int theA;
	G4int theZ;

	// Gamma is A_f(2S_f+1) factor, where A_f is fragment atomic 
	// number and S_f is fragment spin
	G4double Gamma;

	// Discrete Excitation Energies 
	G4RWTValVector<G4double> * ExcitationEnergies;

	//
	G4RWTValVector<G4int> * ExcitationSpins;

};


#endif

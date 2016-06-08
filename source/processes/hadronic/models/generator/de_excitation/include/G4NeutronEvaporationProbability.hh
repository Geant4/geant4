// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999) 
//



#ifndef G4NeutronEvaporationProbability_h
#define G4NeutronEvaporationProbability_h 1


#include "G4EvaporationProbability.hh"


class G4NeutronEvaporationProbability : public G4EvaporationProbability
{
public:
	// Only available constructor
	G4NeutronEvaporationProbability();
		
	~G4NeutronEvaporationProbability() {}
private:  
	// Copy constructor
	G4NeutronEvaporationProbability(const G4NeutronEvaporationProbability &right);

	const G4NeutronEvaporationProbability & operator=(const G4NeutronEvaporationProbability &right);
	G4bool operator==(const G4NeutronEvaporationProbability &right) const;
	G4bool operator!=(const G4NeutronEvaporationProbability &right) const;
  

private:

	virtual G4double CalcAlphaParam(const G4Fragment & fragment) const 
	{ return 0.76+2.2/pow(G4double(fragment.GetA()-GetA()),1.0/3.0);}
	
	virtual G4double CalcBetaParam(const G4Fragment & fragment) const 
	{ return (2.12/pow(G4double(fragment.GetA()-GetA()),2.0/3.0) - 0.05)*MeV/CalcAlphaParam(fragment); }

	// Excitation energy levels 
	G4RWTValVector<G4double> ExcitEnergies;
	// Spin of excitation energy levels 
	G4RWTValVector<G4int> ExcitSpins;

};


#endif

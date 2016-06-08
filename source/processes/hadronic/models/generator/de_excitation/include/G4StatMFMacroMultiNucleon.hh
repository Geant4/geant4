#ifndef G4StatMFMacroMultiNucleon_h
#define G4StatMFMacroMultiNucleon_h 1

#include "G4VStatMFMacroCluster.hh"


class G4StatMFMacroMultiNucleon : public G4VStatMFMacroCluster {

public:

	// Constructor
	G4StatMFMacroMultiNucleon(const G4int Size) : G4VStatMFMacroCluster(Size) {};

	// Destructor
	~G4StatMFMacroMultiNucleon() {};
	

private:

	// Default constructor
	G4StatMFMacroMultiNucleon();
	
	// Copy constructor
	G4StatMFMacroMultiNucleon(const G4StatMFMacroMultiNucleon & right);

	// operators
	G4StatMFMacroMultiNucleon & operator=(const G4StatMFMacroMultiNucleon & right);
	G4bool operator==(const G4StatMFMacroMultiNucleon & right) const;
	G4bool operator!=(const G4StatMFMacroMultiNucleon & right) const;

public:

	G4double CalcMeanMultiplicity(const G4double FreeVol, const G4double mu, 
											const G4double nu, const G4double T);
								
	G4double CalcZARatio(const G4double nu); 
	
	G4double CalcEnergy(const G4double T);

	G4double CalcEntropy(const G4double T, const G4double FreeVol);
	
};

#endif

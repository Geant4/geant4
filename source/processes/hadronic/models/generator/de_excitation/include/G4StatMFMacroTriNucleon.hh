#ifndef G4StatMFMacroTriNucleon_h
#define G4StatMFMacroTriNucleon_h 1

#include "G4VStatMFMacroCluster.hh"
#include "G4NucleiPropertiesTable.hh"

class G4StatMFMacroTriNucleon : public G4VStatMFMacroCluster {

public:

	// Default constructor
	G4StatMFMacroTriNucleon() : G4VStatMFMacroCluster(3) {};

	// Destructor
	~G4StatMFMacroTriNucleon() {};
	

private:

	// Copy constructor
	G4StatMFMacroTriNucleon(const G4StatMFMacroTriNucleon & right);

	// operators
	G4StatMFMacroTriNucleon & operator=(const G4StatMFMacroTriNucleon & right);
	G4bool operator==(const G4StatMFMacroTriNucleon & right) const;
	G4bool operator!=(const G4StatMFMacroTriNucleon & right) const;

public:

	G4double CalcMeanMultiplicity(const G4double FreeVol, const G4double mu, 
											const G4double nu, const G4double T);

	G4double CalcZARatio(const G4double nu) {return theZARatio = 0.5;}

	G4double CalcEnergy(const G4double T);
	
	G4double CalcEntropy(const G4double T, const G4double FreeVol);
	
	
};

#endif

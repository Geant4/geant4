#ifndef G4StatMFMacroBiBiNucleon_h
#define G4StatMFMacroBiNucleon_h 1

#include "G4VStatMFMacroCluster.hh"
#include "G4NucleiPropertiesTable.hh"

class G4StatMFMacroBiNucleon : public G4VStatMFMacroCluster {

public:

	// Default constructor
	G4StatMFMacroBiNucleon() : G4VStatMFMacroCluster(2) {};

	// Destructor
	~G4StatMFMacroBiNucleon() {};
	

private:

	// Copy constructor
	G4StatMFMacroBiNucleon(const G4StatMFMacroBiNucleon & right);

	// operators
	G4StatMFMacroBiNucleon & operator=(const G4StatMFMacroBiNucleon & right);
	G4bool operator==(const G4StatMFMacroBiNucleon & right) const;
	G4bool operator!=(const G4StatMFMacroBiNucleon & right) const;

public:

	G4double CalcMeanMultiplicity(const G4double FreeVol, const G4double mu, 
											const G4double nu, const G4double T);
					
					
					
	G4double CalcZARatio(const G4double nu) {return theZARatio = 0.5;}
	
	G4double CalcEnergy(const G4double T);
	
	G4double CalcEntropy(const G4double T, const G4double FreeVol);
};

#endif

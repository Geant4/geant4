#ifndef G4StatMFMacroTetraNucleon_h
#define G4StatMFMacroTetraNucleon_h 1

#include "G4VStatMFMacroCluster.hh"
#include "G4NucleiPropertiesTable.hh"

class G4StatMFMacroTetraNucleon : public G4VStatMFMacroCluster {

public:

	// Default constructor
	G4StatMFMacroTetraNucleon() : G4VStatMFMacroCluster(4) {};

	// Destructor
	~G4StatMFMacroTetraNucleon() {};
	

private:

	// Copy constructor
	G4StatMFMacroTetraNucleon(const G4StatMFMacroTetraNucleon & right);

	// operators
	G4StatMFMacroTetraNucleon & operator=(const G4StatMFMacroTetraNucleon & right);
	G4bool operator==(const G4StatMFMacroTetraNucleon & right) const;
	G4bool operator!=(const G4StatMFMacroTetraNucleon & right) const;

public:

	G4double CalcMeanMultiplicity(const G4double FreeVol, const G4double mu, 
											const G4double nu, const G4double T);
								
	G4double CalcZARatio(const G4double nu) {return theZARatio = 0.5;}								

	G4double CalcEnergy(const G4double T);

	G4double CalcEntropy(const G4double T, const G4double FreeVol);

};

#endif

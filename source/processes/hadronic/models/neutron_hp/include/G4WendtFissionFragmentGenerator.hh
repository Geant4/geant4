/*
 * File:   G4WendtFissionFragmentGenerator.hh
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on June 21, 2013, 13:58 MST
 */

#include <map>

#include "G4HadFinalState.hh"

#include "G4FissionFragmentGenerator.hh"

class G4WendtFissionFragmentGenerator
{
public:
	G4HadFinalState* ApplyYourself(const G4HadProjectile& projectile, G4int Z, G4int A);
	static G4WendtFissionFragmentGenerator* GetInstance();
	void InitializeANucleus(const G4int A, const G4int Z, const G4int M, const G4String& dataDirectory);
	~G4WendtFissionFragmentGenerator();

private:
	// SINGLETON!!!
	G4WendtFissionFragmentGenerator();
	G4WendtFissionFragmentGenerator(G4WendtFissionFragmentGenerator const&);
	void operator=(G4WendtFissionFragmentGenerator const&);
	// SINGLETON!!!

	/** A map of all the fission isotopes loaded at initialization */
    std::map< const G4int, G4FissionFragmentGenerator* > fissionIsotopes;
    G4NeutronHPNames fileNames;

    G4int Verbosity_;

};

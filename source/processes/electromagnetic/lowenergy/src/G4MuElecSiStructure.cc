//
// G4MuElecSiStructure.cc, 2011/08/29 A.Valentin, M. Raine
//
// Based on the following publications
//
//          - Inelastic cross-sections of low energy electrons in silicon
//	    for the simulation of heavy ion tracks with theGeant4-DNA toolkit,
//	    NSS Conf. Record 2010, p80-85
//	    - Geant4 physics processes for microdosimetry simulation:
//	    very low energy electromagnetic models for electrons in Si,
//	    to be published in TNS
//	    - Geant4 physics processes for microdosimetry simulation:
//	    very low energy electromagnetic models for protons and
//	    heavy ions in Si, to be published in NIMB
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "G4MuElecSiStructure.hh"

G4MuElecSiStructure::G4MuElecSiStructure(): nLevels(6)
{
  energyConstant.push_back(16.65*eV);
  energyConstant.push_back(6.52*eV); 
  energyConstant.push_back(13.63*eV);
  energyConstant.push_back(107.98*eV); 
  energyConstant.push_back(151.55*eV); 
  energyConstant.push_back(1828.5*eV);

  nLevels = energyConstant.size();
}


G4MuElecSiStructure::~G4MuElecSiStructure()
{ }
 

G4double G4MuElecSiStructure::Energy(G4int level)
{
  G4double energ = 0.;

  if (level >=0 && level < nLevels) energ = energyConstant[level];

  return energ;
}

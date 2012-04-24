//
// G4MuElecSiStructure.hh, 2011/08/29 A.Valentin, M. Raine
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


#ifndef G4MUELECSISTRUCTURE_HH
#define G4MUELECSISTRUCTURE_HH 1
 
#include "globals.hh"
#include <vector>

 
class G4MuElecSiStructure
{
public:
  
  G4MuElecSiStructure();
  
  virtual ~G4MuElecSiStructure();
  
  G4double Energy(G4int level);

  G4int NumberOfLevels() { return nLevels; }
  
    
private:
   
 // Number of levels of silicon
 G4int nLevels;

  std::vector<G4double> energyConstant;
  
};

#endif

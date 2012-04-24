//
// G4MuElecElastic.hh, 2011/08/29 A.Valentin, M. Raine
//
// Based on the following publications
//	    - Geant4 physics processes for microdosimetry simulation:
//	    very low energy electromagnetic models for electrons in Si,
//	    to be published in TNS
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef G4MuElecElastic_h
#define G4MuElecElastic_h 1

#include "G4VEmProcess.hh"
#include "G4Electron.hh"

// Available models
#include "G4MuElecElasticModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4MuElecElastic : public G4VEmProcess

{
public: 

  G4MuElecElastic(const G4String& processName ="MuElecElastic",
		     G4ProcessType type = fElectromagnetic);

  virtual ~G4MuElecElastic();

  virtual G4bool IsApplicable(const G4ParticleDefinition&);
  
  virtual void PrintInfo();

protected:

  virtual void InitialiseProcess(const G4ParticleDefinition*);

private:
     
  G4bool       isInitialised;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
#endif

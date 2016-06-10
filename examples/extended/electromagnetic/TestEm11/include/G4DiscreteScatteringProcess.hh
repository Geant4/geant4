#ifndef G4DiscreteScatteringProcess_HH
#define G4DiscreteScatteringProcess_HH

#include "G4VEmProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4DiscreteScatteringProcess : public G4VEmProcess{
        
public:
        
  G4DiscreteScatteringProcess(G4int iNumAngles=1);
  
  virtual ~G4DiscreteScatteringProcess();
  
  virtual void InitialiseProcess(const G4ParticleDefinition*);
  
  virtual bool IsApplicable(const G4ParticleDefinition& p);
  // This will need to check whether the particle is something 
  // our particle works for. Right now, protons.
  
  virtual void PrintInfo();
  
private:
     
  bool  fIsInitialised;    
  G4int fNumAngles;
};




#endif



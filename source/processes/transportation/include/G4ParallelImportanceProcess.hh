#ifndef G4ParallelImportanceProcess_hh
#define G4ParallelImportanceProcess_hh G4ParallelImportanceProcess_hh


#include "G4ParallelTransport.hh"
#include "G4ImportancePostStepDoIt.hh"

class G4VImportanceSampler;
class G4Nsplit_Weight;

class G4ParallelImportanceProcess : public G4ParallelTransport {
public:
  G4ParallelImportanceProcess(const G4VImportanceSampler &aImportanceSampler,
			       G4VPGeoDriver &pgeodriver, 
			       G4VParallelStepper &aStepper,
			       const G4String &aName = "ParallelImportanceProcess");  

  G4VParticleChange *PostStepDoIt(const G4Track& ,
				  const G4Step&);
  
  

private:
  G4ParallelImportanceProcess(const G4ParallelImportanceProcess &);
  G4ParallelImportanceProcess &operator=(const G4ParallelImportanceProcess &);
  
  void Error(const G4String &m){
    G4cout << "ERROE: in G4ImportanceProcess:: " <<  m << G4endl;
    exit(1);
  }

  G4ParticleChange *fParticleChange;
  const G4VImportanceSampler &fImportanceSampler;  
  G4ImportancePostStepDoIt fImportancePostStepDoIt;
};


#endif

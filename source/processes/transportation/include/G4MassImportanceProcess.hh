#ifndef G4MassImportanceProcess_hh
#define G4MassImportanceProcess_hh G4MassImportanceProcess_hh

#include "G4VProcess.hh"
#include "G4ImportancePostStepDoIt.hh"

class G4VImportanceAlgorithm;
class G4ImportanceFinder;
class G4VIStore;

class G4MassImportanceProcess : public G4VProcess {
public:
  G4MassImportanceProcess(const G4VImportanceAlgorithm &aImportanceAlgorithm,
			  const G4VIStore &aIstore,
			  const G4String &aName = "MassImportanceProcess");
  ~G4MassImportanceProcess();
  
  virtual G4double 
  PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				       G4double   previousStepSize,
				       G4ForceCondition* condition);
  
  virtual G4VParticleChange * 
  PostStepDoIt(const G4Track& ,
	       const G4Step&);
  
  //  no operation in  AtRestDoIt and  AlongStepDoIt
  G4double 
  AlongStepGetPhysicalInteractionLength(const G4Track&,
					G4double  ,
					G4double  ,
					G4double& ,
					G4GPILSelection*
					){ return -1.0; };
  
  G4double 
  AtRestGetPhysicalInteractionLength(const G4Track& ,
				     G4ForceCondition* 
				     ) { return -1.0; };
  
  //  no operation in  AtRestDoIt and  AlongStepDoIt
  G4VParticleChange* 
  AtRestDoIt(const G4Track& ,
	     const G4Step&
	     ) {return 0;};
  G4VParticleChange* AlongStepDoIt(const G4Track& ,
				   const G4Step& 
				   ) {return 0;};
  
private:
  G4MassImportanceProcess(const G4MassImportanceProcess &);
  G4MassImportanceProcess &operator=(const G4MassImportanceProcess &);

  G4ParticleChange *fParticleChange;
  G4ImportancePostStepDoIt fImportancePostStepDoIt;
  const G4VImportanceAlgorithm &fImportanceAlgorithm;
  G4ImportanceFinder *fImportanceFinder;
};

#endif

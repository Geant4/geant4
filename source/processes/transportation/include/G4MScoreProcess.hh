#ifndef G4MScoreProcess_hh 
#define G4MScoreProcess_hh G4MScoreProcess_hh

#include "G4VProcess.hh"
class G4VPScorer;

class G4MScoreProcess : public G4VProcess {
public:
  G4MScoreProcess(G4VPScorer &aScorer,
		  const G4String &aName = "MScoreProcess");
  virtual ~G4MScoreProcess();

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
  G4MScoreProcess(const G4MScoreProcess &);
  G4MScoreProcess &operator=(const G4MScoreProcess &);

  G4VPScorer &fScorer;  
};

#endif

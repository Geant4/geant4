#ifndef G4BOptnCloning_hh
#define G4BOptnCloning_hh 1

#include "G4VBiasingOperation.hh"
#include "G4ParticleChange.hh"

// -- very like a splitting by 2, but with weights see freely.

class G4BOptnCloning : public G4VBiasingOperation {
public:
  // -- Constructor :
  G4BOptnCloning(G4String name);
  // -- destructor:
  virtual ~G4BOptnCloning();
  
public:
  // -- Methods from G4VBiasingOperation interface:
  // -------------------------------------------
  // -- Unsed:
  virtual const G4VBiasingInteractionLaw* ProvideOccurenceBiasingInteractionLaw( const G4BiasingProcessInterface* ) {return 0;}
  virtual G4VParticleChange*                             ApplyFinalStateBiasing( const G4BiasingProcessInterface*,
										 const G4Track*,
										 const G4Step*  ) {return 0;}
  // -- Used:
  virtual G4double                                  DistanceToApplyOperation( const G4Track*,
									      G4double,
									      G4ForceCondition* condition)
  {
    *condition = NotForced; return 0; // -- acts immediately
  }
  virtual G4VParticleChange*                    GenerateBiasingFinalState( const G4Track*,
									   const G4Step*  );
  
public:
  // -- Additional methods, specific to this class:
  // ----------------------------------------------
  void SetCloneWeights(G4double clone1Weight, G4double clone2Weight) {fClone1W = clone1Weight ; fClone2W = clone2Weight;}

private:
  G4double fClone1W, fClone2W;
  G4ParticleChange fParticleChange;
};

#endif

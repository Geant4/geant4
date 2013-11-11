#ifndef G4BOptnForceFreeFlight_hh
#define G4BOptnForceFreeFlight_hh 1

#include "G4VBiasingOperation.hh"
#include "G4ForceCondition.hh"
class G4ILawForceFreeFlight;

class G4BOptnForceFreeFlight : public G4VBiasingOperation {
public:
  // -- Constructor :
  G4BOptnForceFreeFlight(G4String name);
  // -- destructor:
  virtual ~G4BOptnForceFreeFlight();
  
public:
  // -- Methods from G4VBiasingOperation interface:
  // -------------------------------------------
  // -- Used:
  virtual const G4VBiasingInteractionLaw* ProvideOccurenceBiasingInteractionLaw( const G4BiasingProcessInterface* );
  virtual G4ForceCondition                                ProposeForceCondition( const G4ForceCondition ) {return Forced;}
  virtual G4bool                                        DenyProcessPostStepDoIt( const G4BiasingProcessInterface*, const G4Track*, const G4Step*, G4double& );
  virtual void                                                      AlongMoveBy( const G4BiasingProcessInterface*, const G4Step*, G4double );

  // -- Unused:
  virtual G4VParticleChange*                       ApplyFinalStateBiasing( const G4BiasingProcessInterface*,
									   const G4Track*,
									   const G4Step*  ) {return 0;}
  virtual G4double                               DistanceToApplyOperation( const G4Track*,
									   G4double,
									   G4ForceCondition*) {return DBL_MAX;}
  virtual G4VParticleChange*                    GenerateBiasingFinalState( const G4Track*,
									   const G4Step*  ) {return 0;}


public:
  // -- Additional methods, specific to this class:
  // ----------------------------------------------
  // -- return concrete type of interaction law:
  G4ILawForceFreeFlight* GetForceFreeFlightLaw() {
    return fForceFreeFlightInteractionLaw;
  }
  // -- initialization for weight:
  void ResetInitialTrackWeight(G4double w) {fInitialTrackWeight = w; fCumulatedWeightChange = 1.0;}
  
private:
  G4ILawForceFreeFlight* fForceFreeFlightInteractionLaw;
  G4double fCumulatedWeightChange, fInitialTrackWeight;
};

#endif

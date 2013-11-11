#ifndef G4BOptrForceCollision_hh
#define G4BOptrForceCollision_hh 1

// Mimic a "force collision" a la MCNP. This is for neutral particles only.
// This is to be considered as an experimental implementation at the beta
// release level.
// When the track enters the volume, it is cloned. One copy makes a forced
// free flight up to the volume exit. The other copy makes a forced collision
// inside the volume. Weights are calculated accordingly.
// Particle processes to be forced must be wrapped with the G4BiasingProcessInterface.
// In this experimental implementation, only gammas and neutrons are considered.

#include "G4VBiasingOperator.hh"
class G4BOptnForceFreeFlight;
class G4BOptnForceCommonTruncatedExp;
class G4BOptnCloning;
class G4VProcess;
class G4BiasingProcessInterface;
class G4ParticleDefinition;
#include <vector>
#include <map>
#include "G4ThreeVector.hh"

class G4BOptrForceCollision : public G4VBiasingOperator {
public:
  G4BOptrForceCollision(G4String particleToForce,                    G4String name="ForceCollision");
  G4BOptrForceCollision(const G4ParticleDefinition* particleToForce, G4String name="ForceCollision");
  ~G4BOptrForceCollision();
  
private:
  // -- Mandatory from base class :
  virtual G4VBiasingOperation*  ProposeOccurenceBiasingOperation(const G4Track* track, const G4BiasingProcessInterface* callingProcess);
  virtual G4VBiasingOperation* ProposeFinalStateBiasingOperation(const G4Track*, const G4BiasingProcessInterface*) {return 0;}
  virtual G4VBiasingOperation* ProposeNonPhysicsBiasingOperation(const G4Track* track, const G4BiasingProcessInterface* callingProcess);
  // -- optional methods from base class:
public:
  virtual void StartTracking( const G4Track* track );
  virtual void   ExitBiasing( const G4Track*, const G4BiasingProcessInterface* );

  using G4VBiasingOperator::OperationApplied; // -- informs compiler of multiple OperationApplied symbols
                                              // -- in base class, to avoid warning messages for hidden functions
                                              // -- in case only one function if overloaded.
  // -- single operation applied (whatever operation):
  virtual void OperationApplied( const G4BiasingProcessInterface*            callingProcess, G4BiasingAppliedCase                      biasingCase,
				 G4VBiasingOperation*                      operationApplied, const G4VParticleChange*        particleChangeProduced );

private:
  std::map< const G4BiasingProcessInterface*, G4BOptnForceFreeFlight* > fFreeFlightOperations;
  G4BOptnForceCommonTruncatedExp*  fSharedForceInteractionOperation;
  G4BOptnCloning*                                 fCloningOperation;
  G4double                                      fInitialTrackWeight;
  //G4double                                         fWeightToRestore;
  std::vector < const G4VProcess* >   fProcesses;
  const G4VProcess *fFirstProcess, *fLastProcess;
  G4bool                                  fSetup;
  const G4ParticleDefinition*          fParticle;
  G4VBiasingOperation* fPreviousOperationApplied;
  const G4Track*                   fTrackToForce;
  G4ThreeVector                fPreviousMomentum;


};

#endif

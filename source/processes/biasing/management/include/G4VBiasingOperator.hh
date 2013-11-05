#ifndef G4VBiasingOperator_hh
#define G4VBiasingOperator_hh 1

#include "globals.hh"

class G4VBiasingOperation;
class G4Track;
class G4BiasingProcessInterface;
class G4LogicalVolume;
class G4VParticleChange;
#include <map>
#include <vector>
#include "G4BiasingAppliedCase.hh"

class G4VBiasingOperator {
public:
  // ---------------
  // -- Constructor:
  // ---------------
  G4VBiasingOperator(G4String name);
  virtual ~G4VBiasingOperator();
  
  // ----------------------------------------------
  // -- abstract and user interface to sub-classes:
  // ----------------------------------------------
protected:
  // -- mandatory methods to let the operator tell about biasing operations to be applied:
  // -------------------------------------------------------------------------------------
  // -- These three methods have the same arguments passed : the current G4Track pointer, and the pointer of the
  // -- G4BiasingProcessInterface instance calling this biasing operator. This same biasing operator will be called by each
  // -- of the G4BiasingProcessInterface instances, meaning for example that:
  // --     - if one G4BiasingProcessInterface with no wrapped physics process exits, ProposeNonPhysicsBiasingOperation(...)
  // --       will be called one time at the beginning of the step,
  // --     - if three G4BiasingProcessInterface instances exist, each of these one wrapping a physics process (eg
  // --       conversion, Compton, photo-electric), ProposeOccurenceBiasingOperation(...) will be called three times,
  // --       by each of these instances, at the beginning of the step and ProposeFinalStateBiasingOperation(...) will
  // --       also be called by each of these instances, at the PostStepDoIt level.
  // -- If a null pointer is returned, the analog -unbiased- behavior is adopted.
  // -- non-physics-based biasing:
  // -----------------------------
  // -- [ First operator method called, at the PostStepGetPhysicalInterationLenght(...) level. ]
  virtual G4VBiasingOperation* ProposeNonPhysicsBiasingOperation( const G4Track* track, const G4BiasingProcessInterface* callingProcess ) = 0;
  // -- physics-based biasing:
  // -------------------------
  // -- Method to propose an occurence biasing operation : ie a change of the interaction length distribution. The proposed
  // -- biasing operation will then be asked for its interaction law.
  // -- Note that *** all sanity checks regarding the operation and its interaction law will have to have been performed
  // -- before returning the biasing operation pointer *** as no corrective/aborting actions will be possible beyond this point.
  // -- The informations provided by the G4BiasingProcessInterface calling process (previous occurence operation, previous step length,
  // -- etc.) might be useful for doing this. They will be useful also to decide with continuing with a same operation proposed
  // -- in the previous step, updating the interaction law taking into account the new G4Track state and the previous step size.
  // -- [ Second operator method called, at the PostStepGetPhysicalInterationLenght(...) level. ]
  virtual G4VBiasingOperation*  ProposeOccurenceBiasingOperation( const G4Track* track, const G4BiasingProcessInterface* callingProcess ) = 0;
  // -- [ Third operator method called, at the PostStepDoIt(...) level. ]
  virtual G4VBiasingOperation* ProposeFinalStateBiasingOperation( const G4Track* track, const G4BiasingProcessInterface* callingProcess ) = 0;
  
protected:
  // -- optionnal methods for further information passed to the operator:
  // --------------------------------------------------------------------
  // ---- report to operator about the operation applied, the biasingCase value provides what sort of biasing applied:
  virtual void OperationApplied( const G4BiasingProcessInterface* callingProcess, G4BiasingAppliedCase               biasingCase,
				 G4VBiasingOperation*           operationApplied, const G4VParticleChange* particleChangeProduced );
  // ---- same as above, report about the operation applied, for the case an occurence biasing was applied, together or not with a final state biasing.
  // ---- The variable biasingCase tells if the final state is a biased one or not. **But in all cases**, this call happens only
  // ---- for an occurence biaising : ie the occurence weight is applied on top of the particleChangeProduced, which is the particle
  // ---- *before* the weight application for occurence biasing.
  virtual void OperationApplied( const G4BiasingProcessInterface*            callingProcess, G4BiasingAppliedCase                      biasingCase,
				 G4VBiasingOperation*             occurenceOperationApplied, G4double                 weightForOccurenceInteraction,
				 G4VBiasingOperation*            finalStateOperationApplied, const G4VParticleChange*        particleChangeProduced );
protected:
  // ---- method to inform operator that its biasing control is over (exit volume, or end of tracking):
  // ---- [Called at the beginning of next step, or at the end of tracking.]
  virtual void ExitBiasing( const G4Track* track, const G4BiasingProcessInterface* callingProcess );


protected:
  // -- Utility:
  // -----------
  // ---- A utility method allowing to store secondaries produced by an operation.
  // ---- The stored secondaries are the ones in the particle change
  void RememberSecondaries( const G4BiasingProcessInterface*         callingProcess,
			    const G4VBiasingOperation*             operationApplied,
			    const G4VParticleChange*         particleChangeProduced );
  // ---- Informations about track is erased:
  void ForgetTrack( const G4Track* track );
  
public:
  // ---- inform the operator of the start of the run:
  virtual void      StartRun() {}
  // ---- inform the operator of the start (end) of the tracking of a new track:
  virtual void StartTracking( const G4Track* /* track */ ) {}
  virtual void   EndTracking() {}
  
  
  
  // --------------------
  // -- public interface:
  // --------------------
  // -- needed by user:
public:
  const G4String                                    GetName() const {return fName;}
  void                                             AttachTo( const G4LogicalVolume* ); // -- attach to single volume 
  void                                         AttachToTree( const G4LogicalVolume* ); // -- attach to mother and full descending tree
  static G4VBiasingOperator*             GetBiasingOperator( const G4LogicalVolume* ); // -- might go to a manager ; or moved to volume
  G4BiasingAppliedCase        GetPreviousBiasingAppliedCase() const {return fPreviousBiasingAppliedCase;}
  // -- all operators (might got to a manager):
  static const std::vector < G4VBiasingOperator* >& GetBiasingOperators() {return fOperators;}
  
  
  // -- used by biasing process interface, or used by an other operator (not expected to be invoked differently than with these two cases):
public:
  G4VBiasingOperation*  GetProposedOccurenceBiasingOperation( const G4Track* track, const G4BiasingProcessInterface* callingProcess );
  G4VBiasingOperation* GetProposedFinalStateBiasingOperation( const G4Track* track, const G4BiasingProcessInterface* callingProcess );
  G4VBiasingOperation* GetProposedNonPhysicsBiasingOperation( const G4Track* track, const G4BiasingProcessInterface* callingProcess );
  void                                        ExitingBiasing( const G4Track* track, const G4BiasingProcessInterface* callingProcess );
  
public:
  void ReportOperationApplied( const G4BiasingProcessInterface* callingProcess, G4BiasingAppliedCase biasingCase,
			       G4VBiasingOperation* operationApplied, const G4VParticleChange* particleChangeProduced );
  void ReportOperationApplied( const G4BiasingProcessInterface* callingProcess, G4BiasingAppliedCase biasingCase,
			       G4VBiasingOperation* occurenceOperationApplied,  G4double weightForOccurenceInteraction,
			       G4VBiasingOperation* finalStateOperationApplied, const G4VParticleChange* particleChangeProduced );
  
  
public:
  const G4VBiasingOperation* GetPreviousNonPhysicsAppliedOperation() {return  fPreviousAppliedNonPhysicsBiasingOperation;}
  const G4VBiasingOperation* GetBirthOperation( const G4Track* );
  
  
private:
  const G4String fName;
  // -- shared, non-thread local:
  static std::map< const G4LogicalVolume*, G4VBiasingOperator* > fLogicalToSetupMap;
  // -- thread local:
  static std::vector < G4VBiasingOperator* > fOperators;


  // -- For this operator:
  std::vector< const G4LogicalVolume* >        fRootVolumes;
  std::map   < const G4LogicalVolume*, G4int > fDepthInTree;
  
  // -- current operation:
  G4VBiasingOperation*     fOccurenceBiasingOperation;
  G4VBiasingOperation*    fFinalStateBiasingOperation;
  G4VBiasingOperation*    fNonPhysicsBiasingOperation;
  
  // -- previous operations:
  const G4VBiasingOperation*     fPreviousProposedOccurenceBiasingOperation;
  const G4VBiasingOperation*    fPreviousProposedFinalStateBiasingOperation;
  const G4VBiasingOperation*    fPreviousProposedNonPhysicsBiasingOperation;
  const G4VBiasingOperation*      fPreviousAppliedOccurenceBiasingOperation;
  const G4VBiasingOperation*     fPreviousAppliedFinalStateBiasingOperation;
  const G4VBiasingOperation*     fPreviousAppliedNonPhysicsBiasingOperation;
  G4BiasingAppliedCase                          fPreviousBiasingAppliedCase;
  
};

#endif

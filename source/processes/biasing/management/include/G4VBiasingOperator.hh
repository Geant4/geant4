//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: $
//
// --------------------------------------------------------------------
// GEANT 4 class header file 
//
// Class Description:
//
//     An abstract class to pilot the biasing in a logical volume. This
// class is for *making decisions* on biasing operations to be applied.
// These ones are represented by the G4VBiasingOperation class.
// The volume in which biasing is applied is specified by the
// AttachTo(const G4LogicalVolume *) method. This has to be specified
// at detector construction time in the method ConstructSDandField() of
// G4VUsedDetectorConstruction.
//
//     At tracking time the biasing operator is messaged by each
// G4BiasingProcessInterface object attached to the current track. For
// example, if three physics processes are under biasing, and if an
// additional G4BiasingProcessInterface is present to handle non-physics
// based biasing (splitting, killing), the operator will be messaged by
// these four G4BiasingProcessInterface objects.
//     The idendity of the calling G4BiasingProcessInterface is known
// to the G4VBiasingOperator by passing this process pointer to the
// operator.
//
// ** Mandatory methods: **
//
// Three types of biasing are to be decided by the G4VBiasingOperator:
//
//   1) non-physics-based biasing:
//   -----------------------------
//   Meant for pure killing/splitting/etc. biasing operations, not 
//   associated to a physics process:
//
//   virtual G4VBiasingOperation* ProposeNonPhysicsBiasingOperation( const G4Track* track,
//                                                                   const G4BiasingProcessInterface* callingProcess ) = 0;
//
//   Arguments are the current track, and the G4BiasingProcessInterface
//   pointer making the call to the operator. In this case, this process
//   does not wrap a physics process and 
//                callingProcess->GetWrappedProcess() == 0.
//
//   The G4VBiasingOperation pointer returned is the operation to be
//   applied. Zero can be returned. This operation will limit the
//   step and propose a final state.
//   
//   This method is the first operator method called, it is called at the
//   by the PostStepGetPhysicalInterationLenght(...) method of the
//   G4BiasingProcessInterface.
//
//   2) physics-based biasing:
//   -------------------------
//   Physics-based biasing operations are of two types:
//     - biasing of the physics process occurence interaction law
//     - biasing of the physics process final state production
//
//   a) The biasing of the occurence interaction law is proposed by:
//
//   virtual G4VBiasingOperation*  ProposeOccurenceBiasingOperation( const G4Track* track, 
//                                                                   const G4BiasingProcessInterface* callingProcess ) = 0;
//   The current G4Track pointer and the G4BiasingProcessInterface
//   pointer of the process calling the operator are passed. The
//   G4BiasingProcessInterface process wraps an actual physics process
//   which pointer can be obtained with
//                callingProcess->GetWrappedProcess() .
//
//   The biasing operation returned will be asked for its biasing
//   interaction by the calling process, which will be a const object
//   for the process. All setup and sampling regarding this law should be done
//   in the operator before returning the related operation to the process.
//
//   This method is the second operator one called in a step, it is called by
//   the PostStepGetPhysicalInterationLenght(...) method of the
//   G4BiasingProcessInterface.
//
//   b) The biasing of the physics process final state is proposed by:
//
//   virtual G4VBiasingOperation* ProposeFinalStateBiasingOperation( const G4Track* track,
//                                                                   const G4BiasingProcessInterface* callingProcess ) = 0;
//
//   The operator can propose a biasing operation that will handle the
//   physic process final state biasing. As in previous case a) the
//   G4BiasingProcessInterface process wraps an actual physics process
//   which pointer can be obtained with:
//                callingProcess->GetWrappedProcess() .
//
//   Cases a) and b) are handled independently, and one or two of these
//   biasing types can be provided in the same step.
//
//   This method is the last operator one called in a step, it is called
//   by the PostStepDoIt(...) method of the G4BiasingProcessInterface.
//
//
// ** Optional methods: **
//
//     At the end of the step, the operator is messaged by the G4BiasingProcessInterface
// for operation(s) which have been applied during the step. One of the two following
// methods is called:
//
// virtual void OperationApplied( const G4BiasingProcessInterface*                callingProcess,
//                                G4BiasingAppliedCase                               biasingCase,
// 				  G4VBiasingOperation*                          operationApplied,
//                                const G4VParticleChange*                particleChangeProduced );
// At most a single biasing operation was applied by the process:
//    - a non-physics biasing operation was applied, biasingCase == BAC_NonPhysics ;
//    - physics-based biasing:
//      - the operator requested no biasing operations, and did let the physics
//        process go : biasingCase ==  BAC_None;
//      - a single final state biasing was proposed, with no concomittant occurence:
//        biasingCase ==  BAC_FinalState;
// The operation applied and final state passed to the tracking (particleChangeProduced) are
// passed as information to the operator.
// 
// virtual void OperationApplied( const G4BiasingProcessInterface*                callingProcess,
//                                G4BiasingAppliedCase                               biasingCase,
//				  G4VBiasingOperation*                 occurenceOperationApplied,
//                                G4double                         weightForOccurenceInteraction,
//				  G4VBiasingOperation*                finalStateOperationApplied,
//                                const G4VParticleChange*                particleChangeProduced );
// This method is called in case an occurence biasing operation has been applied during the step.
// The biasingCase value is then the one of the final state biasing, if any : depending on if the
// occurence operation was applied alone and together with a final state operation, the
// biasingCase will take values:
//     - occurence biasing alone : biasingCase == BAC_None ;
//       in which case finalStateOperationApplied == 0;
//     - occurence biasing + final state biasing : biasingCase ==  BAC_FinalState;
// The particleChangeProduced is the one *before* application of the weight for occurence : hence
// either the particle change of the (analog) physics process, or the biased final state, resulting
// from the biasing by the finalStateOperationApplied operation.
//
//   
//      ----------------G4VBiasingOperation ----------------
//
// Author: M.Verderi (LLR), November 2013
//
// --------------------------------------------------------------------

#ifndef G4VBiasingOperator_hh
#define G4VBiasingOperator_hh 1

#include "globals.hh"

class G4VBiasingOperation;
class G4Track;
class G4BiasingProcessInterface;
class G4LogicalVolume;
class G4VParticleChange;
class G4BiasingOperatorStateNotifier;
#include <map>
#include <vector>
#include "G4BiasingAppliedCase.hh"
#include "G4Cache.hh"


class G4VBiasingOperator {

  // -- State machine used to inform operators
  // -- about run starting.
  // -- Defined at the end of this file.
  friend class G4BiasingOperatorStateNotifier;
  
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
  // -- optional methods for further information passed to the operator:
  // -------------------------------------------------------------------
  // ---- report to operator about the operation applied, the biasingCase value provides the case of biasing applied:
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
  // -----------------------------------
  // -- Delegation to an other operator:
  // -----------------------------------
  // -- An operator may wish to select a sequence of operations already implemented in an
  // -- existing biasing operator. In this case, this operator can delegate its work to
  // -- the "delegated" one by calling DelegateTo( G4VBiasingOperation* delegated );
  // -- §§ Should we have:
  // -- §§    - a "step delegation" -where the delegation is made for the current step only-
  // -- §§    - a long delegation where the delegation can hold over several steps, as long as
  // -- §§      the scheme is not completed. [let's call it "scheme delegation"]
  // -- §§      In this case the "execution/delegated" operator might switch off back the
  // -- §§      delegation from the "delegator" when it knows it has done its work.
  // -- §§ Add a private SetDelegator( G4VBiasingOperator* ) method, call on the delegated
  // -- §§ operator.
  // -- §§ For a step long delegation, the ReportOperationApplied should be used to "unset"
  // -- §§ the delegation. For a scheme long delegation, the delegater operator will unset
  // -- §§ itself has delegation. Likely to happen in the ReportOperationApplied as well,
  // -- §§ but not sure it is mandatory though.


public:
  // ---- Configure() is called in sequential mode or for master thread in MT mode.
  // ---- It is in particular aimed at registering ID's to physics model at run initialization.
  virtual void            Configure() {}
  // ---- ConfigureForWorker() is called in MT mode only, and only for worker threads.
  // ---- It is not not to be used to register ID's to physics model catalog.
  virtual void   ConfigureForWorker() {}
  // ---- inform the operator of the start of the run:
  virtual void              StartRun() {}
  // ---- inform the operator of the start (end) of the tracking of a new track:
  virtual void         StartTracking( const G4Track* /* track */ ) {}
  virtual void           EndTracking() {}
  
  
  
  // --------------------
  // -- public interface:
  // --------------------
  // -- needed by user:
public:
  const G4String                                    GetName() const {return fName;}
  void                                             AttachTo( const G4LogicalVolume* ); // -- attach to single volume 

  G4BiasingAppliedCase        GetPreviousBiasingAppliedCase() const {return fPreviousBiasingAppliedCase;}
  // -- all operators (might got to a manager):
  static const std::vector < G4VBiasingOperator* >& GetBiasingOperators() {return fOperators.Get();}
  // -- get operator associated to a logical volume:
  static G4VBiasingOperator*             GetBiasingOperator( const G4LogicalVolume* ); // -- might go to a manager ; or moved to volume

  
  
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
  
  
private:
  const G4String fName;
  // -- thread local:
  //  static std::map< const G4LogicalVolume*, G4VBiasingOperator* > fLogicalToSetupMap;
  static G4MapCache< const G4LogicalVolume*, G4VBiasingOperator* > fLogicalToSetupMap;
  // -- thread local:
  static G4VectorCache<G4VBiasingOperator* > fOperators;
  // static std::vector < G4VBiasingOperator* > fOperators;

  // -- thread local:
  static G4Cache< G4BiasingOperatorStateNotifier* > fStateNotifier;


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

// -- state machine to get biasing operators
// -- messaged at the beginning of runs:
#include "G4VStateDependent.hh"
class G4BiasingOperatorStateNotifier : public G4VStateDependent {
public:
  G4BiasingOperatorStateNotifier();
  ~G4BiasingOperatorStateNotifier();
public:
  G4bool Notify(G4ApplicationState requestedState);
private:
  G4ApplicationState fPreviousState;
};

#endif

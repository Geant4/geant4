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
//--------------------------------------------------------------------
//
// G4BiasingProcessInterface
//
// Class Description:
//        A wrapper process making the interface between the tracking
//    and the G4VBiasingOperator objects attached to volumes.
//        If this process holds a physics process, it forwards
//    tracking calls to this process in volume where not biasing
//    occurs. In volumes with biasing (with a G4VBiasingOperator
//    attached) the process gets what to do messaging the biasing
//    operator :
//        - at the PostStepGPIL level, for getting an occurence biasing
//          operation. If such an operation is returned to the process
//          this operation will be messaged at several places.
//        - at the PostStepDoIt level, to get a possible final state
//          biasing operation
//        If the process does not hold a physics process, it is meant
//    as handling "non physics" biasing operations: pure splitting
//    or pure kiling for example (ie not brem splitting).
//
//--------------------------------------------------------------------
//   Initial version                         Sep. 2013 M. Verderi

#ifndef G4BiasingProcessInterface_h
#define G4BiasingProcessInterface_h

#include "globals.hh"
#include "G4VProcess.hh"
#include "G4Cache.hh"

class G4VBiasingInteractionLaw;
class G4InteractionLawPhysical;
class G4VBiasingOperator;
class G4VBiasingOperation;
class G4ParticleChangeForOccurenceBiasing;
class G4ParticleChange;
class G4ParticleChangeForNothing;

class G4BiasingProcessInterface : public G4VProcess {
public:
  // --------------------------------------------------------------------------------
  // -- constructor for dealing with biasing options not affecting physics processes:
  // --------------------------------------------------------------------------------
  G4BiasingProcessInterface( G4String name = "biasWrapper(0)" );
  // ---------------------------------------------------------------
  // -- constructor to transform the behaviour of a physics process:
  // ---------------------------------------------------------------
  // -- wrappedProcess pointer MUST NOT be null.
  G4BiasingProcessInterface(G4VProcess* wrappedProcess,
			    G4bool wrappedIsAtRest, G4bool wrappedIsAlongStep, G4bool wrappedIsPostStep,
			    G4String useThisName = "");
  ~G4BiasingProcessInterface();
  // -- pointer of wrapped physics process:
  G4VProcess* GetWrappedProcess() const {return fWrappedProcess;}
  
  // ---------------------
  // -- Biasing interface:
  // ---------------------
  // -- Helper methods:
  // ------------------
  // -- Current step and previous step biasing operator, if any:
  G4VBiasingOperator*              GetCurrentBiasingOperator()  const {return             fCurrentBiasingOperator; }
  G4VBiasingOperator*             GetPreviousBiasingOperator()  const {return            fPreviousBiasingOperator; }
  // -- current and previous operation:
  G4VBiasingOperation*   GetCurrentOccurenceBiasingOperation() const { return          fOccurenceBiasingOperation; }
  G4VBiasingOperation*  GetPreviousOccurenceBiasingOperation() const { return  fPreviousOccurenceBiasingOperation; }
  G4VBiasingOperation*  GetCurrentFinalStateBiasingOperation() const { return         fFinalStateBiasingOperation; }
  G4VBiasingOperation* GetPreviousFinalStateBiasingOperation() const { return fPreviousFinalStateBiasingOperation; }
  G4VBiasingOperation*  GetCurrentNonPhysicsBiasingOperation() const { return         fNonPhysicsBiasingOperation; }
  G4VBiasingOperation* GetPreviousNonPhysicsBiasingOperation() const { return fPreviousNonPhysicsBiasingOperation; }
  


  // ------------------
  // -- Helper methods:
  // ------------------
  // ---- Tell is this process is first/last in the PostStep GPIL and DoIt lists.
  // ---- If physOnly is true, only wrapper for physics processes are considered,
  // ---- otherwise all G4BiasingProcessInterface processes of this particle are
  // ---- considered.
  // ---- These methods just return the corresponding flag values setup at
  // ---- initialization phase by the next four ones.
  // ---- Will not be updated if processes are activate/unactivated on the fly.
  // ---- Use next methods (less fast) instead in this case.
  G4bool GetIsFirstPostStepGPILInterface(G4bool physOnly = true) const;
  G4bool  GetIsLastPostStepGPILInterface(G4bool physOnly = true) const;
  G4bool GetIsFirstPostStepDoItInterface(G4bool physOnly = true) const;
  G4bool  GetIsLastPostStepDoItInterface(G4bool physOnly = true) const;
  // ---- Determine if the process is first/last in the PostStep GPIL and DoIt lists.
  G4bool    IsFirstPostStepGPILInterface(G4bool physOnly = true) const;
  G4bool     IsLastPostStepGPILInterface(G4bool physOnly = true) const;
  G4bool    IsFirstPostStepDoItInterface(G4bool physOnly = true) const;
  G4bool     IsLastPostStepDoItInterface(G4bool physOnly = true) const;
  // -- Information about wrapped process:
  G4bool GetWrappedProcessIsAtRest() const { return fWrappedProcessIsAtRest; }
  G4bool  GetWrappedProcessIsAlong() const { return  fWrappedProcessIsAlong; }
  G4bool   GetWrappedProcessIsPost() const { return   fWrappedProcessIsPost; }
  
 
  // -- Information methods:
  G4double         GetPreviousStepSize()             const { return fPreviousStepSize;}
  G4double       GetCurrentMinimumStep()             const { return fCurrentMinimumStep;}
  G4double           GetProposedSafety()             const { return fProposedSafety;}
  void               SetProposedSafety(G4double sft)       { fProposedSafety = sft;}
  // -- return the actual PostStep and AlongStep limits returned by the process to the tracking :
  G4double             GetPostStepGPIL()             const { return         fBiasingPostStepGPIL; }
  G4double            GetAlongStepGPIL()             const { return fWrappedProcessAlongStepGPIL; }


  // --------------------------------------------------------------
  // --                  G4VProcess interface                    --
  // --------------------------------------------------------------
public:
  // -- Start/End tracking:
  void StartTracking(G4Track* track);
  void   EndTracking();
  
  // -- PostStep methods:
  virtual G4double           PostStepGetPhysicalInteractionLength(const G4Track&                track,
								  G4double           previousStepSize,
								  G4ForceCondition*         condition);
  virtual G4VParticleChange*                         PostStepDoIt(const G4Track&                track,
								  const G4Step&                  step);
  // -- AlongStep methods:
  virtual G4double           AlongStepGetPhysicalInteractionLength(const G4Track&                track,
								   G4double           previousStepSize,
								   G4double         currentMinimumStep, 
								   G4double&            proposedSafety, 
								   G4GPILSelection*          selection);
  virtual G4VParticleChange*                         AlongStepDoIt(const G4Track&                track,
								   const G4Step&                  step);
  // -- AtRest methods
  virtual G4double           AtRestGetPhysicalInteractionLength(const G4Track&,
								G4ForceCondition*);
  virtual G4VParticleChange*                         AtRestDoIt(const G4Track&,
								const G4Step&);
  
  virtual G4bool         IsApplicable(const G4ParticleDefinition& pd);
  virtual void      BuildPhysicsTable(const G4ParticleDefinition& pd);
  virtual void    PreparePhysicsTable(const G4ParticleDefinition& pd);
  virtual G4bool    StorePhysicsTable(const G4ParticleDefinition* pd,
				      const G4String& s, G4bool f);
  virtual G4bool RetrievePhysicsTable(const G4ParticleDefinition* pd,
				      const G4String& s, G4bool f);
  // --
  virtual void SetProcessManager(const G4ProcessManager*); 
  virtual const G4ProcessManager* GetProcessManager(); 
  // --
  virtual void ResetNumberOfInteractionLengthLeft();
  //  virtual void ClearNumberOfInteractionLengthLeft();
  // --
  //  virtual void  DumpInfo() const;
  
  virtual void          SetMasterProcess(G4VProcess* masterP);
  virtual void   BuildWorkerPhysicsTable(const G4ParticleDefinition& pd);
  virtual void PrepareWorkerPhysicsTable(const G4ParticleDefinition& pd);


private:
  // ---- Setup first/last flags:
  void       SetUpFirstLastFlags();
  void  ResetForUnbiasedTracking();

  G4Track*                                               fCurrentTrack;
  G4double                                           fPreviousStepSize;
  G4double                                         fCurrentMinimumStep;
  G4double                                             fProposedSafety;

  G4VBiasingOperator*                              fCurrentBiasingOperator;
  G4VBiasingOperator*                             fPreviousBiasingOperator;
  G4VBiasingOperation*                          fOccurenceBiasingOperation;
  G4VBiasingOperation*                         fFinalStateBiasingOperation;
  G4VBiasingOperation*                         fNonPhysicsBiasingOperation;
  G4VBiasingOperation*                  fPreviousOccurenceBiasingOperation;
  G4VBiasingOperation*                 fPreviousFinalStateBiasingOperation;
  G4VBiasingOperation*                 fPreviousNonPhysicsBiasingOperation;

  G4bool                            fResetWrappedProcessInteractionLength;

  G4VProcess*                                             fWrappedProcess;
  const G4bool                                     fIsPhysicsBasedBiasing;
  const G4bool                                    fWrappedProcessIsAtRest;
  const G4bool                                     fWrappedProcessIsAlong;
  const G4bool                                      fWrappedProcessIsPost;


  G4double                                    fWrappedProcessPostStepGPIL;
  G4double                                           fBiasingPostStepGPIL;
  G4double                               fWrappedProcessInteractionLength; // -- inverse of analog cross-section
  G4ForceCondition                          fWrappedProcessForceCondition;
  G4ForceCondition                                 fBiasingForceCondition;
  G4double                                   fWrappedProcessAlongStepGPIL;
  G4double                                          fBiasingAlongStepGPIL;
  G4GPILSelection                            fWrappedProcessGPILSelection;
  G4GPILSelection                                   fBiasingGPILSelection;

  const G4VBiasingInteractionLaw*                  fBiasingInteractionLaw;
  const G4VBiasingInteractionLaw*          fPreviousBiasingInteractionLaw;
  G4InteractionLawPhysical*                       fPhysicalInteractionLaw;
  G4ParticleChangeForOccurenceBiasing*    fOccurenceBiasingParticleChange;
  G4ParticleChange*                                       fParticleChange; // -- to be changed with a light version for weight only
  G4ParticleChangeForNothing*                        fDummyParticleChange;
  G4bool                                               fFirstLastFlags[8];
  G4int       IdxFirstLast(G4int firstLast, G4int GPILDoIt, G4int physAll) const
  {
    // -- be careful : all arguments are *assumed* to be 0 or 1. No check
    // -- for that is provided. Should be of pure internal usage.
    return 4*firstLast + 2*GPILDoIt + physAll;
  }
  // -- TO DO : let some common operation be done only once, by the first
  // -- process in GPIL. Nothing done yet, may help to reduce overhead.
  // -- For example each process instance fetch for the current operator,
  // -- which could be done once, and with pointer shared.
  // -- Idea is to have a class containing the data members shared among
  // -- the instances.
  // -- Could make these shared data members shared only per particle type
  // -- (at present, with a single static, any instance of a process shares
  // -- its members with any other instance, regardless of particle type.)
  G4bool                                                    fIamFirstGPIL;


  // -- MUST be **thread local**:
  static G4Cache<G4bool>                                     fResetInteractionLaws;
  static G4Cache<G4bool>                                              fCommonStart;
  static G4Cache<G4bool>                                                fCommonEnd;
  
  // -- Maintain lists of interfaces attached to a same particle (ie
  // -- to a same G4ProcessManager ).
  // -- This is used when biasing operations in different processes need to
  // -- cooperate (like for knowing the total cross-section for example).
  std::vector< G4BiasingProcessInterface* >*             fCoInterfaces;
  // -- thread local:
  static G4MapCache< const G4ProcessManager*, 
		     std::vector< G4BiasingProcessInterface* > > fManagerInterfaceMap;
  const G4ProcessManager* fProcessManager;
  
};

#endif

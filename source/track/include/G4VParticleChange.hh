// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VParticleChange.hh,v 1.3 1999-05-06 11:42:52 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD group
//
// ------------------------------------------------------------
//   Implemented for the new scheme                23 Mar. 1998  H.Kurahige
// 
//  This class is the abstract class for ParticleChange.
//
//  The ParticleChange class ontains the results after invocation 
//  of a physics process. This includes final states of parent
//  particle (momentum, energy, etc) and secondary particles generated 
//  by the interaction.
//  The tracking assumes that all the values of energy and
//  momentum are in global reference system, therefore all the
//  needed Lorentz transformations must have been already Done
//  when filling the data-members of this class.
// 
//
//   This abstract class has following four virtual methods
//     virtual G4Step* UpdateStepForAtRest(G4Step* Step);
//     virtual G4Step* UpdateStepForAlongStep(G4Step* Step);
//     virtual G4Step* UpdateStepForPostStep(G4Step* Step);
//     virtual void Initialize(const G4Track&);
//   The UpdateStep methods return the pointer to the G4Step 
//   after updating the given Step information by using final state 
//   information of the track given by a physics process.    
//   User must add methods to keep the final state information 
//   in his derived class as well as implement UpdateStep methods 
//   which he want to use.
//
//   The Initialize methods is provided to refresh the final 
//   state information and should be called by each process 
//   at the beginning of DoIt.
//   
// ------------------------------------------------------------
//   Implement Event Biasing Scheme   9 Nov.,98 H.Kurashige
//   add CheckIt                    13  Apr.,99 H.Kurashige
//   add accuracy leveles            5  May, 99 H.Kurashige
#ifndef G4VParticleChange_h
#define G4VParticleChange_h 1

#include "globals.hh"
#include "G4ios.hh"

class G4Track;
class G4Step;
class G4VEvtBiasMechanism;

#include "G4TrackFastVector.hh"
#include "G4TrackStatus.hh"
#include "G4SteppingControl.hh"


class G4VParticleChange 
{
  public:
    // default constructor
    G4VParticleChange();
    G4VParticleChange(G4bool useEvtBiasing);

    // destructor
    virtual ~G4VParticleChange();

    // equal/unequal operator
    G4bool operator==(const G4VParticleChange &right) const;
    G4bool operator!=(const G4VParticleChange &right) const;
    // "equal" means that teo objects have the same pointer.

  protected:
    // hide copy constructor and assignment operaor as protected
    G4VParticleChange(const G4VParticleChange &right);
    G4VParticleChange & operator=(const G4VParticleChange &right);
 
  public:
    // --- the following methods are for updating G4Step -----   
    virtual G4Step* UpdateStepForAtRest(G4Step* Step);
    virtual G4Step* UpdateStepForAlongStep(G4Step* Step);
    virtual G4Step* UpdateStepForPostStep(G4Step* Step);
    // Return the pointer to the G4Step after updating the Step information
    // by using final state information of the track given by a physics
    // process    
 
  protected:
    G4Step* UpdateStepInfo(G4Step* Step);
    //  Update the G4Step specific attributes 
    //  (i.e. SteppingControl, LocalEnergyDeposit, and TrueStepLength)


  public:
    virtual void Initialize(const G4Track&);
    // This methods will be called by each process at the beginning of DoIt
    // if necessary.

  protected:
    void InitializeTrueStepLength(const G4Track&);
    void InitializeLocalEnergyDeposit(const G4Track&);
    void InitializeSteppingControl(const G4Track&);
    void InitializeParentWeight(const G4Track&);

    void InitializeStatusChange(const G4Track&);
    void InitializeSecondaries(const G4Track&);
   // ------------------------------------------------------   

  public:
    //---- the following methods are for TruePathLength ----
    G4double GetTrueStepLength() const;
    void  SetTrueStepLength(G4double truePathLength);
    //  Get/Set theTrueStepLength

    //---- the following methods are for LocalEnergyDeposit ----   
    G4double GetLocalEnergyDeposit() const;
    void SetLocalEnergyDeposit(G4double anEnergyPart);
    //  Get/Set the locally deposited energy 

    //---- the following methods are for TrackStatus -----   
    G4TrackStatus GetStatusChange() const;
    void SetStatusChange(G4TrackStatus status); 
    //  Get/Set the final TrackStatus of the current particle.
    // ------------------------------------------------------   

    //---- the following methods are for managements of SteppingControl --
    G4SteppingControl GetSteppingControl() const;
    void SetSteppingControl(G4SteppingControl StepControlFlag);
    //  Set/Get a flag to control stepping manager behavier 
    // ------------------------------------------------------   

    //---- the following methods are for managements of secondaries --
    void Clear();
    //  Clear the contents of this objects 
    //  This method should be called after the Tracking(Stepping) 
    //  manager removes all secondaries in theListOfSecondaries 

    void SetNumberOfSecondaries(G4int totSecondaries);
    //  SetNumberOfSecondaries must be called just before AddSecondary()
    //  in order to secure memory space for theListOfSecondaries 
    //  This method resets theNumberOfSecondaries to 0
    //  (that will be incremented at every AddSecondary() call).

    G4int GetNumberOfSecondaries() const;
    //  Returns the number of secondaries current stored in
    //  G4TrackFastVector.

    G4Track* GetSecondary(G4int anIndex) const;
    //  Returns the pointer to the generated secondary particle
    //  which is specified by an Index.

    void AddSecondary(G4Track* aSecondary);
    //  Add a secondary particle to theListOfSecondaries.
    // ------------------------------------------------------   

    G4double GetParentWeight() const ;

    virtual void DumpInfo() const;
    //  Print out information

    void SetVerboseLevel(G4int vLevel);
    G4int GetVerboseLevel() const;

  protected:

    G4TrackFastVector* theListOfSecondaries;
    //  The vector of secondaries.

    G4int theNumberOfSecondaries;
    //  The total number of secondaries produced by each process.

    G4int theSizeOftheListOfSecondaries;
    //  TheSizeOftheListOfSecondaries;

    G4TrackStatus theStatusChange;
    //  The changed (final) track status of a given particle.

    G4SteppingControl theSteppingControlFlag;     
    //  a flag to control stepping manager behavior 

    G4double theLocalEnergyDeposit;
    //  It represents the part of the energy lost for discrete
    //  or semi-continuous processes which is due to secondaries
    //  not generated because they would have been below their cut
    //  threshold.
    //  The sum of the locally deposited energy + the delta-energy
    //  coming from the continuous processes gives the
    //  total energy loss localized in the current Step.

    G4double theTrueStepLength;
    //  The value of "True" Step Length
    
    G4int verboseLevel;
    //  The Verbose level

  public:
    // CheckIt method is provided for debug
    virtual G4bool CheckIt(const G4Track&);
 
    // CheckIt method is activated 
    // if debug flag is set and 'G4VERBOSE' is defined 
    void   ClearDebugFlag();
    void   SetDebugFlag();
    G4bool GetDebugFlag() const; 
        
  protected: 
    G4bool   debugFlag;
 
    // accuracy levels
    static const G4double accuracyForWarning;
    static const G4double accuracyForException; 

  //---- following methods and members are used for Event Biasing
  public:
    virtual void   RegisterEBMechanism(G4VEvtBiasMechanism* );
    virtual void   SwOnEB();
    virtual void   SwOffEB();
    virtual G4bool IsEBActive() const;
    virtual G4VEvtBiasMechanism* GetEBMechanism();

  protected:
    G4VEvtBiasMechanism* theEBMechanism;
    G4bool fUseEB;
    G4double theParentWeight;
     
};

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VEvtBiasMechanism.hh"
#include "G4VParticleChange.icc"

#endif

















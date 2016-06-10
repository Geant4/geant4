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
// $Id: G4VParticleChange.hh 68795 2013-04-05 13:24:46Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//
// ------------------------------------------------------------
//   Implemented for the new scheme                23 Mar. 1998  H.Kurahige
//
// Class Description 
//  This class is the abstract class for ParticleChange.
//-
//  The ParticleChange class ontains the results after invocation 
//  of a physics process. This includes final states of parent
//  particle (momentum, energy, etc) and secondary particles generated 
//  by the interaction.
//  The tracking assumes that all the values of energy and
//  momentum are in global reference system, therefore all the
//  needed Lorentz transformations must have been already Done
//  when filling the data-members of this class.
//-
//-
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
//-
//   The Initialize methods is provided to refresh the final 
//   state information and should be called by each process 
//   at the beginning of DoIt.
//   
// ------------------------------------------------------------
//   Implement Event Biasing Scheme   9 Nov.,98 H.Kurashige
//   add CheckIt                    13  Apr.,99 H.Kurashige
//   add accuracy leveles            5  May, 99 H.Kurashige
//   add check secondaries          11  June, 03 H.Kurashige
//   add new methods of ProposeXXX  08  May, 04 H.Kurashige  
//   remove obsolete methods of SetXXX  19  Sep, 04 H.Kurashige  
//   add flag for first/last step in volume 30 Oct. 2006 H.Kurashige
//   add nonIonizingEnergyLoss          26 Mar 2007 H.Kurashige 
//   modify/fix bugs related to weight  17 Sep. 2011   H.Kurashige 
//   fix bugs related to weight         29 Apr. 2012   H.Kurashige
//

#ifndef G4VParticleChange_h
#define G4VParticleChange_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <cmath>

class G4Track;
class G4Step;

#include "G4TrackFastVector.hh"
#include "G4TrackStatus.hh"
#include "G4SteppingControl.hh"


class G4VParticleChange 
{
  public:
    // default constructor
    G4VParticleChange();

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
 
  public: // with description
    // --- the following methods are for updating G4Step -----   
    virtual G4Step* UpdateStepForAtRest(G4Step* Step);
    virtual G4Step* UpdateStepForAlongStep(G4Step* Step);
    virtual G4Step* UpdateStepForPostStep(G4Step* Step);
    // Return the pointer to the G4Step after updating the Step information
    // by using final state information of the track given by a physics
    // process    
 
  protected: // with description
    G4Step* UpdateStepInfo(G4Step* Step);
    //  Update the G4Step specific attributes 
    //  (i.e. SteppingControl, LocalEnergyDeposit, and TrueStepLength)


  public: // with description
    virtual void Initialize(const G4Track&);
    // This methods will be called by each process at the beginning of DoIt
    // if necessary.

  protected:
    void InitializeTrueStepLength(const G4Track&);
    void InitializeLocalEnergyDeposit(const G4Track&);
    void InitializeSteppingControl(const G4Track&);
    void InitializeParentWeight(const G4Track&);
    void InitializeParentGlobalTime(const G4Track&);

    void InitializeStatusChange(const G4Track&);
    void InitializeSecondaries(const G4Track&);
    void InitializeStepInVolumeFlags(const G4Track&);
   // ------------------------------------------------------   
 
  public: // with description
    //---- the following methods are for TruePathLength ----
    G4double GetTrueStepLength() const;
    void  ProposeTrueStepLength(G4double truePathLength);
    //  Get/Propose theTrueStepLength

    //---- the following methods are for LocalEnergyDeposit ----   
    G4double GetLocalEnergyDeposit() const;
    void ProposeLocalEnergyDeposit(G4double anEnergyPart);
    //  Get/Propose the locally deposited energy 
 
    //---- the following methods are for nonIonizingEnergyDeposit  ----   
    G4double GetNonIonizingEnergyDeposit() const;
    void ProposeNonIonizingEnergyDeposit(G4double anEnergyPart);
    //  Get/Propose the non-ionizing deposited energy 

    //---- the following methods are for TrackStatus -----   
    G4TrackStatus GetTrackStatus() const;
    void ProposeTrackStatus(G4TrackStatus status); 
    //  Get/Propose the final TrackStatus of the current particle.
    // ------------------------------------------------------   

    //---- the following methods are for managements of SteppingControl --
    G4SteppingControl GetSteppingControl() const;
    void ProposeSteppingControl(G4SteppingControl StepControlFlag);
    //  Set/Propose a flag to control stepping manager behavier 
    // ------------------------------------------------------   
 
    //---- the following methods are for managements of initial/last step
    G4bool GetFirstStepInVolume() const;
    G4bool GetLastStepInVolume() const;
    void   ProposeFirstStepInVolume(G4bool flag);
    void   ProposeLastStepInVolume(G4bool flag);

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

    G4double GetWeight() const;
    G4double GetParentWeight() const ;
    //  Get weight of the parent (i.e. current) track
    void ProposeWeight(G4double finalWeight);
    void ProposeParentWeight(G4double finalWeight);
    //  Propse new weight of the parent (i.e. current) track
    //  As for AlongStepDoIt, the parent weight will be set 
    //  in accumulated manner
    //  i.e.) If two processes propose weight of W1 and W2 respectively
    //  for the track with initial weight of W0 
    //  the final weight is set to
    //  (W1/W0) * (W2/W0) * W0  
  
    void     SetSecondaryWeightByProcess(G4bool);
    G4bool   IsSecondaryWeightSetByProcess() const;  
    // In default (fSecondaryWeightByProcess flag is false), 
    // the weight of secondary tracks will be set to 
    // the parent weight
    // If fSecondaryWeightByProcess flag is true, 
    // the weight of secondary tracks will not be changed 
    // by the ParticleChange
    // (i.e. the process determine the secodary weight)
    // NOTE: 
    // Make sure that only one processe in AlongStepDoIt 
    // proposes the parent weight, 
    // If several processes in AlongStepDoIt proposes 
    // the parent weight and add secondaties with 
    // fSecondaryWeightByProcess is set to false, 
    // secondary weights may be wrong

    void   SetParentWeightByProcess(G4bool);
    G4bool   IsParentWeightSetByProcess() const;  
    // Obsolete

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

    G4double theNonIonizingEnergyDeposit;
    //   non-ionizing energu deposit is defined as 
    //   a part of local energy deposit, which does not cause
    //   ionization of atoms
 
    G4double theTrueStepLength;
    //  The value of "True" Step Length
    

    G4bool theFirstStepInVolume;
    G4bool theLastStepInVolume;
    // flag for initial/last step

    G4double theParentWeight;
    // Weight ofparent track
    G4bool isParentWeightProposed;
    // flags for Weight ofparent track
    G4bool   fSetSecondaryWeightByProcess;  
    //  flag for setting weight of secondaries  
 
    G4double theParentGlobalTime;
    // global time of the parent. 
    // This is used only for checking

    G4int verboseLevel;
    //  The Verbose level

  public: // with description
    // CheckIt method is provided for debug
    virtual G4bool CheckIt(const G4Track&);
 
    // CheckIt method is activated 
    // if debug flag is set and 'G4VERBOSE' is defined 
    void   ClearDebugFlag();
    void   SetDebugFlag();
    G4bool GetDebugFlag() const; 

  protected:
    // CheckSecondary method is provided for debug
    G4bool CheckSecondary(G4Track&);
 
    G4double GetAccuracyForWarning() const;
    G4double GetAccuracyForException() const;

  protected: 
    G4bool   debugFlag;
 
    // accuracy levels
    static const G4double accuracyForWarning;
    static const G4double accuracyForException; 


};

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VParticleChange.icc"

#endif

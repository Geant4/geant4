// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRestContinuousDiscreteProcess.hh,v 1.1 1999-01-07 16:13:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD group
// ------------------------------------------------------------
//   New Physics scheme            8  Mar. 1997  H.Kurahige
// ------------------------------------------------------------
//   modified                     26  Mar. 1997  H.Kurahige
//   modified                     16  Apr. 1997  L.Urban    
//   modified AlongStepGPIL etc.  17 Dec. 1997 H.Kurashige
//   fix bugs in GetGPILSelection() 24 Jan. 1998 H.Kurashige
//   modified for new ParticleChange 12 Mar. 1998  H.Kurashige


#ifndef G4VRestContinuousDiscreteProcess_h
#define G4VRestContinuousDiscreteProcess_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VProcess.hh"

class G4VRestContinuousDiscreteProcess : public G4VProcess 
{
  //  Abstract class which defines the public behavior of
  //  discrete physics interactions.
  public:     

      G4VRestContinuousDiscreteProcess(const G4String& ,
				       G4ProcessType   aType = fNotDefined );
      G4VRestContinuousDiscreteProcess(G4VRestContinuousDiscreteProcess &);

      ~G4VRestContinuousDiscreteProcess();

      G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double   previousStepSize,
			     G4ForceCondition* condition
			    );

      G4VParticleChange* PostStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    );

      G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double previousStepSize,
			     G4double currentMinimumStep,
			     G4double& currentSafety,
                             G4GPILSelection* selection
			    );

      G4VParticleChange* AlongStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    );
 
      G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& ,
			     G4ForceCondition* 
			    );

      G4VParticleChange* AtRestDoIt(
			     const G4Track& ,
			     const G4Step&
			    );

   protected:
      virtual G4double GetMeanLifeTime(const G4Track& aTrack,G4ForceCondition* condition)=0;
      //  Calculates the mean life-time (i.e. for decays) of the
      //  particle at rest due to the occurence of the given process,
      //  or converts the probability of interaction (i.e. for
      //  annihilation) into the life-time of the particle for the
      //  occurence of the given process.

    virtual G4double GetContinuousStepLimit(const G4Track& aTrack,
                             G4double  previousStepSize,
                             G4double  currentMinimumStep,
			     G4double& currentSafety
                                                             )=0;
  private:
    // this is the returnd value of  G4GPILSelection in 
    // the arguments of AlongStepGPIL()
    G4GPILSelection  valueGPILSelection;

  protected:
    // these two methods are set/get methods for valueGPILSelection
    // these two methods are set/get methods for valueGPILSelection
    void SetGPILSelection(G4GPILSelection selection)
    { valueGPILSelection = selection;};

    G4GPILSelection GetGPILSelection() const{return valueGPILSelection;};

  protected:
     virtual G4double GetMeanFreePath(const G4Track& aTrack,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                                                               )=0;
      //  Calculates from the macroscopic cross section a mean
      //  free path, the value is returned in units of distance.

 private:
  // hide default constructor and assignment operator as private 
      G4VRestContinuousDiscreteProcess();
      G4VRestContinuousDiscreteProcess & operator=(const G4VRestContinuousDiscreteProcess &right);

};

// -----------------------------------------
//  inlined function members implementation
// -----------------------------------------
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4MaterialTable.hh"
#include "G4VParticleChange.hh"

inline G4double G4VRestContinuousDiscreteProcess::AtRestGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4ForceCondition* condition
			    )
{
  // beggining of tracking 
  ResetNumberOfInteractionLengthLeft();

  // condition is set to "Not Forced"
  *condition = NotForced;

  // get mean life time
  currentInteractionLength = GetMeanLifeTime(track, condition);

#ifdef G4VERBOSE
   if ((currentInteractionLength <0.0) || (verboseLevel>2)){
    G4cout << "G4VRestContinuousDiscreteProcess::AtRestGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" <<endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " << track.GetMaterial()->GetName() <<endl;
    G4cout << "MeanLifeTime = " << currentInteractionLength/ns << "[ns]" <<endl;
  }
#endif

  return theNumberOfInteractionLengthLeft * currentInteractionLength;
}


inline G4VParticleChange* G4VRestContinuousDiscreteProcess::AtRestDoIt( 
			     const G4Track&,
			     const G4Step& 
			    )
{
//  clear NumberOfInteractionLengthLeft
    ClearNumberOfInteractionLengthLeft();

    return pParticleChange;
}

inline G4double G4VRestContinuousDiscreteProcess::AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double previousStepSize,
			     G4double currentMinimumStep,
			     G4double& currentSafety,
                             G4GPILSelection* selection
			    )
{
  // GPILSelection is set to defaule value of CandidateForSelection
  valueGPILSelection = CandidateForSelection;

  // get Step limit proposed by the process
  G4double steplength = GetContinuousStepLimit(track,previousStepSize,currentMinimumStep, currentSafety);

  // set return value for G4GPILSelection
  *selection = valueGPILSelection;

#ifdef G4VERBOSE
  if (verboseLevel>1){
    G4cout << "G4VRestContinuousDiscreteProcess::AlongStepGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" <<endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " <<  track.GetMaterial()->GetName() <<endl;
    G4cout << "IntractionLength= " << steplength/cm <<"[cm] " <<endl;
  }
#endif
  return  steplength ;
}

inline G4VParticleChange* G4VRestContinuousDiscreteProcess::AlongStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    )
{ 
    return pParticleChange;
}

inline G4double G4VRestContinuousDiscreteProcess::PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double   previousStepSize,
			     G4ForceCondition* condition
			    )
{
  if ( (previousStepSize <=0.0) || (theNumberOfInteractionLengthLeft<=0.0)) {
    // beggining of tracking (or just after DoIt of this process)
    ResetNumberOfInteractionLengthLeft();
  } else {
    // subtract NumberOfInteractionLengthLeft 
    SubtractNumberOfInteractionLengthLeft(previousStepSize);
    if(theNumberOfInteractionLengthLeft<perMillion)
       theNumberOfInteractionLengthLeft=0.;         
  }

  // condition is set to "Not Forced"
  *condition = NotForced;

  // get mean free path
  currentInteractionLength = GetMeanFreePath(track, previousStepSize, condition);

#ifdef G4VERBOSE
   if ((currentInteractionLength <=0.0) || (verboseLevel>2)){
    G4cout << "G4VRestContinuousDiscreteProcess::PostStepGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" <<endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " << track.GetMaterial()->GetName() <<endl;
    G4cout << "MeanFreePath = " << currentInteractionLength/cm << "[cm]" <<endl;
  }
#endif

  G4double value = theNumberOfInteractionLengthLeft * currentInteractionLength;
#ifdef G4VERBOSE
   if (verboseLevel>1){
    G4cout << "G4VRestContinuousDiscreteProcess::PostStepGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" <<endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " <<  track.GetMaterial()->GetName() <<endl;
    G4cout << "InteractionLength= " << value/cm <<"[cm] " <<endl;
  }
#endif
  return value;
}

inline G4VParticleChange* G4VRestContinuousDiscreteProcess::PostStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    )
{ 
//  clear NumberOfInteractionLengthLeft
    ClearNumberOfInteractionLengthLeft();

    return pParticleChange;
}


#endif


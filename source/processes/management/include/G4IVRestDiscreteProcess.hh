// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IVRestDiscreteProcess.hh,v 1.6 1999-12-15 14:53:42 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $Id: 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD group
//
// Class Description  
//  Abstract class which defines the public behavior of
//  rest + discrete physics interactions.
//
// ------------------------------------------------------------
//   New Physics scheme           8  Mar. 1997  H.Kurahige
// ------------------------------------------------------------


#ifndef G4IVRestDiscreteProcess_h
#define G4IVRestDiscreteProcess_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VProcess.hh"


class G4IVRestDiscreteProcess : public G4VProcess 
{
  //  Abstract class which defines the public behavior of
  //  rest + discrete physics interactions.

  public: // with description    

    // constructors 
     G4IVRestDiscreteProcess(const G4String& ,
			    G4ProcessType   aType = fNotDefined );
     G4IVRestDiscreteProcess(G4IVRestDiscreteProcess &);

  public:
     virtual ~G4IVRestDiscreteProcess();

  public: // with description    
     // GPIL and DoIt methods derived from G4VProcess
     virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double   previousStepSize,
			     G4ForceCondition* condition
			    );

     virtual G4VParticleChange* PostStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    );

     virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& ,
			     G4ForceCondition* 
			    );
      
     G4VParticleChange* AtRestDoIt(
			     const G4Track& ,
			     const G4Step&
			    );

     //  no operation in  AlongStepDoIt
     virtual G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track&,
			     G4double  ,
			     G4double  ,
			     G4double& ,
                             G4GPILSelection*
			    ){ return -1.0; }

     //  no operation in  AlongStepDoIt
     virtual G4VParticleChange* AlongStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    ) {return 0;}
 
  protected:// with description    
     virtual void SubtractNumberOfInteractionLengthLeft(
                             G4double previousStepSize) ;

     virtual G4double GetMeanLifeTime(const G4Track& aTrack,G4ForceCondition* condition)=0;
      //  Calculates the mean life-time (i.e. for decays) of the
      //  particle at rest due to the occurence of the given process,
      //  or converts the probability of interaction (i.e. for
      //  annihilation) into the life-time of the particle for the
      //  occurence of the given process.

  private:
  // hide default constructor and assignment operator as private 
      G4IVRestDiscreteProcess();
      G4IVRestDiscreteProcess & operator=(const G4IVRestDiscreteProcess &right);


   protected:
      G4PhysicsTable* theNlambdaTable ;
      G4PhysicsTable* theInverseNlambdaTable ;

      const G4double BIGSTEP;
};

// -----------------------------------------
//  inlined function members implementation
// -----------------------------------------
#include "Randomize.hh"              
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4MaterialTable.hh"
#include "G4VParticleChange.hh"

inline 
 void     G4IVRestDiscreteProcess::SubtractNumberOfInteractionLengthLeft(
                             G4double )
 {
  // dummy routine
   ;
 }    

inline G4VParticleChange* G4IVRestDiscreteProcess::PostStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    )
{ 
//  reset NumberOfInteractionLengthLeft
    ClearNumberOfInteractionLengthLeft();

    return pParticleChange;
}

inline G4double G4IVRestDiscreteProcess::AtRestGetPhysicalInteractionLength(
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
    G4cout << "G4IVRestDiscreteProcess::AtRestGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" <<G4endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " << track.GetMaterial()->GetName() <<G4endl;
    G4cout << "MeanLifeTime = " << currentInteractionLength/ns << "[ns]" <<G4endl;
  }
#endif

  return theNumberOfInteractionLengthLeft * currentInteractionLength;
}


inline G4VParticleChange* G4IVRestDiscreteProcess::AtRestDoIt( 
			     const G4Track&,
			     const G4Step& 
			    )
{
//  clear NumberOfInteractionLengthLeft
    ClearNumberOfInteractionLengthLeft();

    return pParticleChange;
}



#endif














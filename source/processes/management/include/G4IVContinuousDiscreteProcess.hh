// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IVContinuousDiscreteProcess.hh,v 1.3 1999-04-17 06:24:57 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $Id: 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD group
// ------------------------------------------------------------
//   New Physics scheme            8  Mar. 1997  H.Kurahige
// ------------------------------------------------------------
//   modified                     26 Mar. 1997 H.Kurashige
//   modified                     16 Apr. 1997 L.Urban     
//   modified AlongStepGPIL etc.  17 Dec. 1997 H.Kurashige
//   fix bugs in GetGPILSelection() 24 Jan. 1998 H.Kurashige
//   modified for new ParticleChange 12 Mar. 1998  H.Kurashige

#ifndef G4IVContinuousDiscreteProcess_h
#define G4IVContinuousDiscreteProcess_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VProcess.hh"

class G4IVContinuousDiscreteProcess : public G4VProcess 
{
  //  Abstract class which defines the public behavior of
  //  discrete physics interactions.
  public:     

      G4IVContinuousDiscreteProcess(const G4String& ,
				   G4ProcessType   aType = fNotDefined );
      G4IVContinuousDiscreteProcess(G4IVContinuousDiscreteProcess &);

      virtual ~G4IVContinuousDiscreteProcess();

      virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double   previousStepSize,
			     G4ForceCondition* condition
			    );

      virtual G4VParticleChange* PostStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    );

      virtual G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track&,
                             G4double  previousStepSize,
                             G4double  currentMinimumStep,
			     G4double& currentSafety,
                             G4GPILSelection* selection
			    )  ;

      virtual G4VParticleChange* AlongStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    );
 
     //  no operation in  AtRestDoIt
      virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& ,
			     G4ForceCondition* 
			    ) { return -1.0; };

     //  no operation in  AtRestDoIt
      virtual G4VParticleChange* AtRestDoIt(
			     const G4Track& ,
			     const G4Step&
			    ) {return 0;};

  protected:
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
    //------------------------------------------------------
    virtual void SubtractNumberOfInteractionLengthLeft(
                                   G4double previousStepSize) ;  

    // these two methods are set/get methods for valueGPILSelection
    void SetGPILSelection(G4GPILSelection selection)
    { valueGPILSelection = selection;};

    G4GPILSelection GetGPILSelection() const{return valueGPILSelection;};

   private:
  // hide default constructor and assignment operator as private 
      G4IVContinuousDiscreteProcess();
      G4IVContinuousDiscreteProcess & operator=(const G4IVContinuousDiscreteProcess &right);

   protected:
      G4PhysicsTable* theNlambdaTable ; 
      G4PhysicsTable* theInverseNlambdaTable ; 
      const G4double BIGSTEP ;
};
// -----------------------------------------
//  inlined function members implementation
// -----------------------------------------
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4EnergyLossTables.hh"
#include "G4MaterialTable.hh"


inline void G4IVContinuousDiscreteProcess::
                             SubtractNumberOfInteractionLengthLeft(
                             G4double 
                            )
{
 // dummy routine
  ;
}  


inline G4VParticleChange* G4IVContinuousDiscreteProcess::PostStepDoIt(
			     const G4Track& ,
			     const G4Step&
			    )
{ 
//  clear  NumberOfInteractionLengthLeft
    ClearNumberOfInteractionLengthLeft();
    return pParticleChange;
}

inline G4VParticleChange* G4IVContinuousDiscreteProcess::AlongStepDoIt(
			     const G4Track& ,
			     const G4Step&
			    )
{ 
//  clear  NumberOfInteractionLengthLeft
    ClearNumberOfInteractionLengthLeft();
    return pParticleChange;
}

inline G4double G4IVContinuousDiscreteProcess::AlongStepGetPhysicalInteractionLength(
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

  if (verboseLevel>1){
    G4cout << "G4IVContinuousDiscreteProcess::AlongStepGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" <<endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " <<  track.GetMaterial()->GetName() <<endl;
    G4cout << "IntractionLength= " << steplength/cm <<"[cm] " <<endl;
  }
  return  steplength ;
}

#endif








// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VDiscreteProcess.hh,v 1.3 1999-11-07 17:11:46 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD group
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
//	add G4VDiscreteProcess(const G4String&) 24 Jul 1996, Hisaya kurashige
//
// Class Description  
//  Abstract class which defines the public behavior of
//  discrete physics interactions.
//
// ------------------------------------------------------------
//   New Physics scheme           18 Dec. 1996  H.Kurahige
// ------------------------------------------------------------
//   modified                     25 Feb. 1997  H.Kurahige
//   modified                     8  Mar. 1997  H.Kurahige
//   modified                     26 Mar. 1997  H.Kurahige
//   modified                     16 Apr. 1997  L.Urban    
//   modified                     17 Dec. 1997  L.Urban    
//   modified for new ParticleChange 12 Mar. 1998  H.Kurashige

#ifndef G4VDiscreteProcess_h
#define G4VDiscreteProcess_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VProcess.hh"

class G4VDiscreteProcess : public G4VProcess 
{
  //  Abstract class which defines the public behavior of
  //  discrete physics interactions.
  public:     

      G4VDiscreteProcess(const G4String& ,
			 G4ProcessType   aType = fNotDefined );
      G4VDiscreteProcess(G4VDiscreteProcess &);

      virtual ~G4VDiscreteProcess();

  public :// with description
      virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double   previousStepSize,
			     G4ForceCondition* condition
			    );

      virtual G4VParticleChange* PostStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    );

     //  no operation in  AtRestDoIt and  AlongStepDoIt
      virtual G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track&,
			     G4double  ,
			     G4double  ,
			     G4double& ,
                             G4GPILSelection*
			    ){ return -1.0; };

      virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& ,
			     G4ForceCondition* 
			    ) { return -1.0; };

     //  no operation in  AtRestDoIt and  AlongStepDoIt
      virtual G4VParticleChange* AtRestDoIt(
			     const G4Track& ,
			     const G4Step&
			    ) {return 0;};

      virtual G4VParticleChange* AlongStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    ) {return 0;};
 
  protected:// with description
     virtual G4double GetMeanFreePath(const G4Track& aTrack,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                                                               )=0;
      //  Calculates from the macroscopic cross section a mean
      //  free path, the value is returned in units of distance.
 
  private:
  // hide default constructor and assignment operator as private 
      G4VDiscreteProcess();
      G4VDiscreteProcess & operator=(const G4VDiscreteProcess &right);

};

// -----------------------------------------
//  inlined function members implementation
// -----------------------------------------
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4MaterialTable.hh"
#include "G4VParticleChange.hh"

inline G4double G4VDiscreteProcess::PostStepGetPhysicalInteractionLength(
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
    G4cout << "G4VDiscreteProcess::PostStepGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" <<endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " << track.GetMaterial()->GetName() <<endl;
    G4cout << "MeanFreePath = " << currentInteractionLength/cm << "[cm]" <<endl;
  }
#endif


  G4double value = theNumberOfInteractionLengthLeft * currentInteractionLength;
#ifdef G4VERBOSE
  if (verboseLevel>1){
    G4cout << "G4VDiscreteProcess::PostStepGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" <<endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " <<  track.GetMaterial()->GetName() <<endl;
    G4cout << "InteractionLength= " << value/cm <<"[cm] " <<endl;
  }
#endif
  return value;
}

inline G4VParticleChange* G4VDiscreteProcess::PostStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    )
{ 
//  clear NumberOfInteractionLengthLeft
    ClearNumberOfInteractionLengthLeft();

    return pParticleChange;
}


#endif




// This code implementation is the intellectual property of
// the  GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRestProcess.hh,v 1.3 1999-11-07 17:11:48 kurasige Exp $
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
// 
// Class Description 
//  Abstract class which defines the public behavior of
//  physics interactions at rest.
//
// ------------------------------------------------------------
//   New Physics scheme           18 Dec. 1996  H.Kurahige
// ------------------------------------------------------------
//   modified                     25 Feb. 1997  H.Kurahige
//   modified                      8 Mar. 1997  H.Kurahige
//   modified                     26 Mar. 1997  H.Kurahige
//   modified                     16 Apr. 1997  L.Urban    
//   modified                     17 Dec. 1997  L.Urban    
//   modified for new ParticleChange 12 Mar. 1998  H.Kurashige


#ifndef G4VRestProcess_h
#define G4VRestProcess_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VProcess.hh"


class G4VRestProcess : public G4VProcess 
{
  //  Abstract class which defines the public behavior of
  //  physics interactions at rest.

  public:
      G4VRestProcess(const G4String&  ,
		     G4ProcessType   aType = fNotDefined );
      G4VRestProcess(G4VRestProcess& );

      virtual ~G4VRestProcess();

  public:   //  with description
      virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4ForceCondition* condition
			    );

      virtual G4VParticleChange* AtRestDoIt(
			     const G4Track& ,
			     const G4Step&
			    );

     //  no operation in  PostStepDoIt and  AlongStepDoIt
      virtual G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track&,
			     G4double  ,
			     G4double  ,
			     G4double& ,
	                     G4GPILSelection*
			   ){ return -1.0; };

      virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& ,
			     G4double   ,
			     G4ForceCondition* 
			    ) { return -1.0; };

     //  no operation in  PostStepDoIt and  AlongStepDoIt
      virtual G4VParticleChange* PostStepDoIt(
			     const G4Track& ,
			     const G4Step&
			    ) {return 0;};

      virtual G4VParticleChange* AlongStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    ) {return 0;};
 
  protected: //  with description

      virtual G4double GetMeanLifeTime(const G4Track& aTrack,G4ForceCondition* condition)=0;
      //  Calculates the mean life-time (i.e. for decays) of the
      //  particle at rest due to the occurence of the given process,
      //  or converts the probability of interaction (i.e. for
      //  annihilation) into the life-time of the particle for the
      //  occurence of the given process.

 private:
  // hide default constructor and assignment operator as private 
      G4VRestProcess();
      G4VRestProcess & operator=(const G4VRestProcess &right);
};

// -----------------------------------------
//  inlined function members implementation
// -----------------------------------------
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4MaterialTable.hh"
#include "G4VParticleChange.hh"

inline G4double G4VRestProcess::AtRestGetPhysicalInteractionLength(
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
    G4cout << "G4VRestProcess::AtRestGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" <<endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " << track.GetMaterial()->GetName() <<endl;
    G4cout << "MeanLifeTime = " << currentInteractionLength/ns << "[ns]" <<endl;
  }
#endif 

  return theNumberOfInteractionLengthLeft * currentInteractionLength;
}


inline G4VParticleChange* G4VRestProcess::AtRestDoIt( 
			     const G4Track&,
			     const G4Step& 
			    )
{
//  clear NumberOfInteractionLengthLeft
    ClearNumberOfInteractionLengthLeft();

    return pParticleChange;
}


#endif




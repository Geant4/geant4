// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IVRestDiscreteProcess.hh,v 1.2 1999-04-09 10:36:18 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $Id: 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD group
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
  public:     

     G4IVRestDiscreteProcess(const G4String& ,
			    G4ProcessType   aType = fNotDefined );
     G4IVRestDiscreteProcess(G4IVRestDiscreteProcess &);

     ~G4IVRestDiscreteProcess();

     G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double   previousStepSize,
			     G4ForceCondition* condition
			    );

     G4VParticleChange* PostStepDoIt(
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

     //  no operation in  AlongStepDoIt
     G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track&,
			     G4double  ,
			     G4double  ,
			     G4double& ,
                             G4GPILSelection*
			    ){ return -1.0; };

     //  no operation in  AlongStepDoIt
     G4VParticleChange* AlongStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    ) {return NULL;};
 
  protected:
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


};

// -----------------------------------------
//  inlined function members implementation
// -----------------------------------------
#include "Randomize.hh"              
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4MaterialTable.hh"
#include "G4VParticleChange.hh"
#include "G4EnergyLossTables.hh"

inline 
 void     G4IVRestDiscreteProcess::SubtractNumberOfInteractionLengthLeft(
                             G4double )
 {
  // dummy routine
   ;
 }    


inline G4double G4IVRestDiscreteProcess::
                             PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            )
{// get particle,particle type,kin.energy,material,mat.index
  G4double nl,nlold,range,rangeold,rangenext,
           KineticEnergyOld,KineticEnergyNext,value;
  G4bool isOut;
  const G4double BIGSTEP=1.e10 ;
  const G4DynamicParticle* particle = track.GetDynamicParticle();
  const G4ParticleDefinition* particletype = particle->GetDefinition() ;
  G4double KineticEnergy = particle->GetKineticEnergy();
  G4Material* material = track.GetMaterial();
  const G4MaterialTable* theMaterialTable =
                         G4Material::GetMaterialTable();
  G4int materialindex = material->GetIndex();

  nl = (*theNlambdaTable)[materialindex]->
                           GetValue(KineticEnergy,isOut);
  range = G4EnergyLossTables::GetPreciseRangeFromEnergy(particletype,
                                       KineticEnergy,material) ;

  if ( (previousStepSize <=0.0) || (theNumberOfInteractionLengthLeft<=0.0)) {
    // beggining of tracking (or just after DoIt of this process)
    ResetNumberOfInteractionLengthLeft();
  } else {
    // subtract NumberOfInteractionLengthLeft

    rangeold = range + previousStepSize ;
    KineticEnergyOld = G4EnergyLossTables::GetPreciseEnergyFromRange(
                                             particletype,
                                       rangeold,material);
    nlold = (*theNlambdaTable)[materialindex]->
                           GetValue(KineticEnergyOld,isOut);

    if(nlold < nl) {
#ifdef G4VERBOSE
      if(verboseLevel>2) {
	G4cout << GetProcessName() << " PostStepGPIL : Nlambda has been" <<
	  " increased at update.Nlambda old/new :" << nlold <<
	    "  " << nl << endl;
	G4cout << "(theNumberOfInteractionLengthLeft has been increased!)" << endl;
	
	G4cout << " correction : Nlambda old=new ........." << endl;
	
      }
#endif
      //corr. of num errror
      nlold = nl ;
     }

    theNumberOfInteractionLengthLeft -= nlold-nl ;

    if(theNumberOfInteractionLengthLeft<perMillion)
       theNumberOfInteractionLengthLeft=0.;
  }


 // condition is set to "Not Forced"
  *condition = NotForced;

  if(nl <= theNumberOfInteractionLengthLeft){
    value = BIGSTEP ;
  } else {
    KineticEnergyNext = (*theInverseNlambdaTable)[materialindex]->
                        GetValue(nl-theNumberOfInteractionLengthLeft,isOut);
    rangenext = G4EnergyLossTables::GetPreciseRangeFromEnergy(particletype,
                                       KineticEnergyNext,material);

    value = range - rangenext ;

    if(range<rangenext) {
#ifdef G4VERBOSE
      if(verboseLevel>2) {
	G4cout << GetProcessName() << " PostStepGPIL: Step < 0.!, Step=" <<
	         value << endl;
	G4cout << "range,rangenext:" << range << "  " << rangenext << endl ;
	G4cout << "correction : rangenext=range ....." << endl;
      }
#endif 
      //corr. of num error
      rangenext = range ;
      value = range - rangenext ;
    }
    
  }

  return value;
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
    G4cout << "[ " << GetProcessName() << "]" <<endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " << track.GetMaterial()->GetName() <<endl;
    G4cout << "MeanLifeTime = " << currentInteractionLength/ns << "[ns]" <<endl;
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


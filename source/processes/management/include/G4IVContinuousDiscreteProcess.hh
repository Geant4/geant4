// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IVContinuousDiscreteProcess.hh,v 1.1 1999-01-07 16:13:52 gunter Exp $
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

      ~G4IVContinuousDiscreteProcess();

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
                             const G4Track&,
                             G4double  previousStepSize,
                             G4double  currentMinimumStep,
			     G4double& currentSafety,
                             G4GPILSelection* selection
			    )  ;

      G4VParticleChange* AlongStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    );
 
     //  no operation in  AtRestDoIt
      G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& ,
			     G4ForceCondition* 
			    ) { return -1.0; };

     //  no operation in  AtRestDoIt
      G4VParticleChange* AtRestDoIt(
			     const G4Track& ,
			     const G4Step&
			    ) {return NULL;};

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
      G4double BIGSTEP ;
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

inline G4double G4IVContinuousDiscreteProcess::
                             PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            )
{// get particle,particle type,kin.energy,material,mat.index
  G4double nl,nlold,range,rangeold,rangenext,
           KineticEnergyOld,KineticEnergyNext,value;
  G4bool isOut;
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







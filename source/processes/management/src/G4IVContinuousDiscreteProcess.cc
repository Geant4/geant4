// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IVContinuousDiscreteProcess.cc,v 1.2 1999-04-17 06:14:59 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $Id: 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD Group
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
// --------------------------------------------------------------
//   New Physics scheme           8 Jan. 1997  H.Kurahige
// ------------------------------------------------------------

#include "G4IVContinuousDiscreteProcess.hh"
G4IVContinuousDiscreteProcess::G4IVContinuousDiscreteProcess()
                   :G4VProcess("No Name Discrete Process") 
{
  G4Exception("G4IVContinuousDiscreteProcess:: default constructor is called");
}

G4IVContinuousDiscreteProcess::G4IVContinuousDiscreteProcess(const G4String& aName , G4ProcessType aType)
                  : G4VProcess(aName, aType)
{
}

G4IVContinuousDiscreteProcess::~G4IVContinuousDiscreteProcess()
{
}

G4IVContinuousDiscreteProcess::G4IVContinuousDiscreteProcess(G4IVContinuousDiscreteProcess& right)
                  : G4VProcess(right)
{
}


G4double G4IVContinuousDiscreteProcess::
                              PostStepGetPhysicalInteractionLength(
                              const G4Track& track,
                              G4double   previousStepSize,
                              G4ForceCondition* condition
                             )
 {
// get particle,particle type,kin.energy,material,mat.index
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
      G4cout << GetProcessName() <<" PostStepGPIL: Step < 0.!, Step=" << value <<endl;
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








// Implementation of class for selecting ionisation model depending on
// logical volume
//
//

#include "G4IonisationByLogicalVolume.hh"

///////////////////////////////////////////////////////////////////////////
//
//

G4IonisationByLogicalVolume::
G4IonisationByLogicalVolume( const G4String& particleName,
                                   G4LogicalVolume* volumeForPAImodel,
                             const G4String& processName  ) :
  G4VContinuousDiscreteProcess(processName),
  fPAIonisation(NULL),feIonisation(NULL),fMuIonisation(NULL),fhIonisation(NULL)
{
  fParticleName = particleName ;
  fVolumeForPAImodel =  volumeForPAImodel ;
  fMaterialNameForPAI = volumeForPAImodel->GetMaterial()->GetName() ;

  fPAIonisation = new G4PAIonisation(fMaterialNameForPAI);

  if ( fParticleName == "e+" || fParticleName == "e-" )
  {
    feIonisation = new G4eIonisation() ;
  }
  else if ( fParticleName == "mu+" || fParticleName == "mu-" )
  {
    fMuIonisation = new G4MuIonisation() ;
  }
  else
  {
    fhIonisation = new G4hIonisation() ;
  }

}

///////////////////////////////////////////////////////////////////////////
//
//

G4IonisationByLogicalVolume::~G4IonisationByLogicalVolume() 
{
  if(fPAIonisation) delete fPAIonisation ;
  if(feIonisation)  delete feIonisation  ;
  if(fMuIonisation) delete fMuIonisation ;
  if(fhIonisation)  delete fhIonisation  ;

}

///////////////////////////////////////////////////////////////////////////
//
// Methods

G4bool G4IonisationByLogicalVolume::
IsApplicable( const G4ParticleDefinition& particle )
{
  return( particle.GetPDGCharge() != 0.) ;
}

/////////////////////////////////////////////////////////////////////////
//
//

/* ***********************************************
  
G4double G4IonisationByLogicalVolume::
GetConstraints(const G4DynamicParticle* aParticle,
                     G4Material*        aMaterial              )
{
  if ( aMaterial->GetName() == fMaterialNameForPAI )
  {
    return fPAIonisation->GetConstraints(aParticle,aMaterial) ;
  }
  else
  {
    if ( fParticleName == "e+" || fParticleName == "e-" )
    {
      return feIonisation->GetConstraints(aParticle,aMaterial) ;
    }
    else if ( fParticleName == "mu+" || fParticleName == "mu-" )
    {
      return fMuIonisation->GetConstraints(aParticle,aMaterial) ;
    }
    else
    {
      return fhIonisation->GetConstraints(aParticle,aMaterial) ;
    }
  }
}

************************************** */

////////////////////////////////////////////////////////////////////////////
//
//
 
G4VParticleChange* 
G4IonisationByLogicalVolume::PostStepDoIt( const G4Track& track,         
                                           const G4Step&  step      ) 
{
  if ( track.GetVolume()->GetLogicalVolume() == fVolumeForPAImodel )
  {
    return fPAIonisation->PostStepDoIt(track,step) ;
  }
  else
  {
    if ( fParticleName == "e+" || fParticleName == "e-" )
    {
      return feIonisation->PostStepDoIt(track,step) ;
    }
    else if ( fParticleName == "mu+" || fParticleName == "mu-" )
    {
      return fMuIonisation->PostStepDoIt(track,step) ;
    }
    else
    {
      return fhIonisation->PostStepDoIt(track,step) ;
    }
  }
}


////////////////////////////////////////////////////////////////////////
//
//


void G4IonisationByLogicalVolume::
BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
{
  fPAIonisation->BuildPhysicsTable(aParticleType) ;

  if ( fParticleName == "e+" || fParticleName == "e-" )
  {
    return feIonisation->BuildPhysicsTable(aParticleType) ;
  }
  else if ( fParticleName == "mu+" || fParticleName == "mu-" )
  {
    return fMuIonisation->BuildPhysicsTable(aParticleType) ;
  }
  else
  {
    return fhIonisation->BuildPhysicsTable(aParticleType) ;
  }
}

////////////////////////////////////////////////////////////////////////
//
//

G4double G4IonisationByLogicalVolume::
GetContinuousStepLimit( const G4Track& track,
                              G4double previousStepSize,
                              G4double currentMinimumStep,
                              G4double& currentSafety      )  
{
  if ( track.GetVolume()->GetLogicalVolume() == fVolumeForPAImodel )
  {
    return fPAIonisation->GetContinuousStepLimit(track,previousStepSize,
                                      currentMinimumStep,currentSafety) ;
  }
  else
  {
    if ( fParticleName == "e+" || fParticleName == "e-" )
    {
      return feIonisation->GetContinuousStepLimit(track,previousStepSize,
                                      currentMinimumStep,currentSafety) ;
    }
    else if ( fParticleName == "mu+" || fParticleName == "mu-" )
    {
      return fMuIonisation->GetContinuousStepLimit(track,previousStepSize,
                                      currentMinimumStep,currentSafety) ;
    }
    else
    {
      return fhIonisation->GetContinuousStepLimit(track,previousStepSize,
                                      currentMinimumStep,currentSafety) ;
    }
  }
}

////////////////////////////////////////////////////////////////////////
//
//

G4double G4IonisationByLogicalVolume::
GetMeanFreePath( const G4Track&          track,
                       G4double          previousStepSize,
                       G4ForceCondition* condition          ) 
{
  if ( track.GetVolume()->GetLogicalVolume() == fVolumeForPAImodel )
  {
    return fPAIonisation->GetMeanFreePath(track,previousStepSize,condition) ;
  }
  else
  {
    if ( fParticleName == "e+" || fParticleName == "e-" )
    {
      return feIonisation->GetMeanFreePath(track,previousStepSize,condition) ;
    }
    else if ( fParticleName == "mu+" || fParticleName == "mu-" )
    {
      return fMuIonisation->GetMeanFreePath(track,previousStepSize,condition) ;
    }
    else
    {
      return fhIonisation->GetMeanFreePath(track,previousStepSize,condition) ;
    }
  }
}

/////////////////////////////////////////////////////////////////////////
//
//

G4VParticleChange* 
G4IonisationByLogicalVolume::AlongStepDoIt( const G4Track& track ,
                                            const G4Step&  step      ) 
{
  if ( track.GetVolume()->GetLogicalVolume() == fVolumeForPAImodel )
  {
    return fPAIonisation->AlongStepDoIt(track,step) ;
  }
  else
  {
    if ( fParticleName == "e+" || fParticleName == "e-" )
    {
      return feIonisation->AlongStepDoIt(track,step) ;
    }
    else if ( fParticleName == "mu+" || fParticleName == "mu-" )
    {
      return fMuIonisation->AlongStepDoIt(track,step) ;
    }
    else
    {
      return fhIonisation->AlongStepDoIt(track,step) ;
    }
  }
}

//////////////////////////////////////////////////////////////////////////
//
//

void G4IonisationByLogicalVolume::
SetVolumeForPAImodel( G4LogicalVolume* volumeForPAImodel  )
{
  fVolumeForPAImodel = volumeForPAImodel ;
}



//
//
/////////////////////////////////////////////////////////////////////////


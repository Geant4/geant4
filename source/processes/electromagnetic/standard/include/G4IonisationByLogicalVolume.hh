// It is process which allows to select PAIonisation model for some logic volume
// in an application. In other logical volumes standard ionisation models
// depending on particle type are realised.
// The process is declared in PhysicsList::ConstructEM() instead of declaration
// of ionisation processes
//
//
//
// History:
// 04.01.00 V. Grichine, first version based on recommendations of J.Apostolakis
//                       and S. Giani

#ifndef G4IonisationByLogicalVolume_h
#define G4IonisationByLogicalVolume_h 1
 
#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4VPAIenergyLoss.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Electron.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsVector.hh"
#include "G4LogicalVolume.hh"

#include "G4VContinuousDiscreteProcess.hh"
#include "G4PAIonisation.hh"
#include "G4eIonisation.hh"
#include "G4MuIonisation.hh"
#include "G4hIonisation.hh"


class G4IonisationByLogicalVolume : public G4VContinuousDiscreteProcess
{
public:
       
  G4IonisationByLogicalVolume( const G4String& particleName,
                                     G4LogicalVolume* volumeForPAImodel,
                               const G4String& processName );

  ~G4IonisationByLogicalVolume() ;


  // Methods

  G4bool IsApplicable(const G4ParticleDefinition&);
  
  //  G4double GetConstraints(const G4DynamicParticle *aParticle,
  //                                G4Material *aMaterial              );
 
  G4VParticleChange* PostStepDoIt( const G4Track& track,         
                                   const G4Step& Step      ) ;                 

  void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);

  G4double 
  GetContinuousStepLimit( const G4Track& track,
                                G4double previousStepSize,
                                G4double currentMinimumStep,
                                G4double& currentSafety      ) ; 

  G4double GetMeanFreePath( const G4Track& track,
                                  G4double previousStepSize,
                                  G4ForceCondition* condition ) ;

  G4VParticleChange* AlongStepDoIt( const G4Track& track ,
                                    const G4Step& Step      ) ;



   void SetVolumeForPAImodel( G4LogicalVolume* volumeForPAImodel  );



protected:

private:

  G4String fParticleName ;
  G4String fMaterialNameForPAI ;
  G4LogicalVolume* fVolumeForPAImodel ;

  G4bool fTriggerPAI ;

  G4PAIonisation*  fPAIonisation ;
  G4eIonisation*   feIonisation ;
  G4MuIonisation*  fMuIonisation ;
  G4hIonisation*   fhIonisation ;


} ;

#endif

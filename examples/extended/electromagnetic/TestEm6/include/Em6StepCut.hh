// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Class Description:
// Step cut is defined
// Class Description - end
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em6StepCut_h
#define Em6StepCut_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4Step.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em6StepCut : public G4VDiscreteProcess
{
public: // Without description   

     Em6StepCut(const G4String& processName ="UserStepCut" );
     Em6StepCut(Em6StepCut &);

    ~Em6StepCut();

     G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double   previousStepSize,
			     G4ForceCondition* condition
			    );

     G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

     void SetMaxStep(G4double);

  protected:

     // it is not needed here !
     G4double GetMeanFreePath(const G4Track& aTrack,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            );

			    
  private:
  
    // hide assignment operator as private 
    Em6StepCut& operator=(const Em6StepCut &right);

  private:

     G4double MaxChargedStep;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// inlined function members implementation

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4EnergyLossTables.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double Em6StepCut::PostStepGetPhysicalInteractionLength(
                             const G4Track& aTrack,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            )
{
  // condition is set to "Not Forced"
  *condition = NotForced;
   G4double ProposedStep = DBL_MAX;

   if((MaxChargedStep > 0.) &&
      (aTrack.GetVolume() != NULL) &&
      (aTrack.GetVolume()->GetName() == "Absorber") &&
      (aTrack.GetDynamicParticle()->GetDefinition()->GetPDGCharge() != 0.))
        ProposedStep = MaxChargedStep;

   return ProposedStep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VParticleChange* Em6StepCut::PostStepDoIt(
                             const G4Track& aTrack,
                             const G4Step&
                            )
{
   // do nothing
   aParticleChange.Initialize(aTrack);
   return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double Em6StepCut::GetMeanFreePath(const G4Track& aTrack,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            )
{
  return 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif


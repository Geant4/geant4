// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em5StepCut.hh,v 1.3 2001-01-09 12:36:37 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em5StepCut_h
#define Em5StepCut_h 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4Step.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em5StepCut : public G4VDiscreteProcess
{
  public:     

     Em5StepCut(const G4String& processName ="UserStepCut" )
         : G4VDiscreteProcess(processName),MaxChargedStep(DBL_MAX) {};

    ~Em5StepCut() {};
    
     void SetMaxStep(G4double step) {MaxChargedStep = step;};
     
     G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double   previousStepSize,
			     G4ForceCondition* condition
			    );

     G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

     G4double GetMeanFreePath(const G4Track& aTrack,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            ) {return 0.;};     // it is not needed here !

  private:

     G4double MaxChargedStep;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//
// inlined function members implementation

inline G4double Em5StepCut::PostStepGetPhysicalInteractionLength(
                             const G4Track& aTrack,
                             G4double   previousStepSize,
                             G4ForceCondition* condition )
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

inline G4VParticleChange* Em5StepCut::PostStepDoIt(const G4Track& aTrack,
                                                   const G4Step&        )
{
   // do nothing
   aParticleChange.Initialize(aTrack);
   return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif


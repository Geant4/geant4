//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#ifndef hTestStepCut_h
#define hTestStepCut_h 1

//---------------------------------------------------------------------------
//
//---------------------------------------------------------------------------
//
// ClassName:   hTestStepCut
//  
// Description: Cut on step length of charged particles
//
// Authors:   08.04.01  V.Ivanchenko 
//
// Modified:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4Step.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestStepCut : public G4VDiscreteProcess
{
public: // Without description   

   hTestStepCut(const G4String& processName = "UserStepCut" );
   hTestStepCut(hTestStepCut &);

  ~hTestStepCut();

   G4double PostStepGetPhysicalInteractionLength(const G4Track&, G4double,
			                         G4ForceCondition*);

   G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

   inline void SetMaxStep(G4double val) {maxChargedStep = val;};

protected:

     // it is not needed here !
   G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);
			    
private:
  
  // hide assignment operator as private 
  hTestStepCut(const hTestStepCut &right);
  const hTestStepCut& operator=(const hTestStepCut &right);

  G4double maxChargedStep;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// inlined function members implementation

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4EnergyLossTables.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double hTestStepCut::PostStepGetPhysicalInteractionLength(
                        const G4Track& aTrack,
                              G4double   previousStepSize,
                              G4ForceCondition* condition)
{
  // condition is set to "Not Forced"
  *condition = NotForced;
   G4double step = DBL_MAX;

   if((maxChargedStep > 0.) &&
      (aTrack.GetVolume() != 0) &&
      (aTrack.GetVolume()->GetName() == "Absorber") &&
      (aTrack.GetDynamicParticle()->GetDefinition()->GetPDGCharge() != 0.))
        step = maxChargedStep;

   return step;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VParticleChange* hTestStepCut::PostStepDoIt(const G4Track& aTrack,
                                                     const G4Step&)
{
   // do nothing
   aParticleChange.Initialize(aTrack);
   return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double hTestStepCut::GetMeanFreePath(const G4Track&, G4double,
                                              G4ForceCondition*)
{
  return 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif


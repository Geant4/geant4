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
#ifndef test31StepCut_h
#define test31StepCut_h 1

//---------------------------------------------------------------------------
//
//---------------------------------------------------------------------------
//
// ClassName:   test31StepCut
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
#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "test31PhysicsList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class test31PhysicsList;

class test31StepCut : public G4VDiscreteProcess
{
public: // Without description   

   test31StepCut(const G4String&, const test31PhysicsList*);

  ~test31StepCut();

   G4double PostStepGetPhysicalInteractionLength(const G4Track&, G4double,
			                         G4ForceCondition*);

   G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

protected:

     // it is not needed here !
  G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);
			    
private:
  
  // hide assignment operator as private 
  test31StepCut(const test31StepCut &right);
  const test31StepCut& operator=(const test31StepCut &right);

  const test31PhysicsList* thePhysics;
};

// inlined function members implementation

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double test31StepCut::PostStepGetPhysicalInteractionLength(
                        const G4Track& aTrack,
                              G4double,
                              G4ForceCondition* condition)
{
  // condition is set to "Not Forced"
  *condition = NotForced;
   G4double step = DBL_MAX;

   if(aTrack.GetVolume()->GetName() == "Absorber")
        step = thePhysics->GetMaxChargedStep();

   return step;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VParticleChange* test31StepCut::PostStepDoIt(const G4Track& aTrack,
                                                     const G4Step&)
{
   // do nothing
   aParticleChange.Initialize(aTrack);
   return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double test31StepCut::GetMeanFreePath(const G4Track&, G4double,
                                              G4ForceCondition*)
{
  return 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif


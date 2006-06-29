//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
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
#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "hTestPhysicsList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestPhysicsList;

class hTestStepCut : public G4VDiscreteProcess
{
public: // Without description   

   hTestStepCut(const G4String&, const hTestPhysicsList*);

  ~hTestStepCut();

   G4double PostStepGetPhysicalInteractionLength(const G4Track&, G4double,
			                         G4ForceCondition*);

   G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

protected:

     // it is not needed here !
  G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);
			    
private:
  
  // hide assignment operator as private 
  hTestStepCut(const hTestStepCut &right);
  const hTestStepCut& operator=(const hTestStepCut &right);

  const hTestPhysicsList* thePhysics;
};

// inlined function members implementation

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double hTestStepCut::PostStepGetPhysicalInteractionLength(
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


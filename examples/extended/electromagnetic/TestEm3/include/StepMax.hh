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
// $Id: StepMax.hh,v 1.6 2005/03/01 17:55:19 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef StepMax_h
#define StepMax_h 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"

#include "DetectorConstruction.hh"

class StepMaxMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class StepMax : public G4VDiscreteProcess
{
public:

  StepMax(const G4String& processName = "UserStepMax");
 ~StepMax();

  G4bool   IsApplicable(const G4ParticleDefinition&);

  void     SetStepMax(G4int, G4double);

  G4double GetStepMax(G4int k) {return stepMax[k];};

  G4double PostStepGetPhysicalInteractionLength( const G4Track& track,
			                       G4double previousStepSize,
			                       G4ForceCondition* condition);

  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  G4double GetMeanFreePath(const G4Track&, G4double,G4ForceCondition*)
     {return DBL_MAX;};    

private:

  G4double stepMax[MaxAbsor];
     
  StepMaxMessenger* pMess;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline 
G4double StepMax::PostStepGetPhysicalInteractionLength( const G4Track& aTrack,
                                    G4double, G4ForceCondition* condition)
{
  // condition is set to "Not Forced"
  *condition = NotForced;
  
  G4double limit = DBL_MAX; 
  G4int n = aTrack.GetVolume()->GetCopyNo();
  if (n < MaxAbsor) limit = stepMax[n];
  return limit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline 
G4VParticleChange* StepMax::PostStepDoIt(const G4Track& aTrack, const G4Step&)
{
   // do nothing
   aParticleChange.Initialize(aTrack);
   return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


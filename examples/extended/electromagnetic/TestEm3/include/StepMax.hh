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
/// \file electromagnetic/TestEm3/include/StepMax.hh
/// \brief Definition of the StepMax class
//
// $Id: StepMax.hh 98762 2016-08-09 14:08:07Z gcosmo $
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

  G4double GetStepMax(G4int k) { return fStepMax[k];};

  G4double PostStepGetPhysicalInteractionLength( const G4Track& track,
                                               G4double previousStepSize,
                                               G4ForceCondition* condition);

  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  G4double GetMeanFreePath(const G4Track&, G4double,G4ForceCondition*)
     {return DBL_MAX;};    

private:

  G4double fStepMax[kMaxAbsor];
     
  StepMaxMessenger* fMess;
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
  if (n < kMaxAbsor) limit = fStepMax[n];
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


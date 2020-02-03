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
/// \file electromagnetic/TestEm11/include/StepMax.hh
/// \brief Definition of the StepMax class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef StepMax_h
#define StepMax_h 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"

class StepMaxMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class StepMax : public G4VDiscreteProcess
{
  public:

     StepMax(const G4String& processName = "UserMaxStep",
                   G4ProcessType type    = fUserDefined);
    ~StepMax();

     virtual G4bool IsApplicable(const G4ParticleDefinition&);

     void SetMaxStep1(G4double);
     void SetMaxStep2(G4double);     
     void ApplyMaxStep2(G4bool);

     virtual G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                                               G4double previousStepSize,
                                               G4ForceCondition* condition);

     virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

     virtual G4double GetMeanFreePath(const G4Track&,G4double,G4ForceCondition*)
     {return DBL_MAX;};    

  private:

     G4double fMaxStep1;
     G4double fMaxStep2;
     G4bool   fApplyMaxStep2;
     
     StepMaxMessenger* fMess;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


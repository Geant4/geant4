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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "Run.hh"

#include "G4RunManager.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Proton.hh"
#include "G4DNAGenericIonsManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction()
:G4UserSteppingAction()
{
  fSumOfStepLength    =  0;
  fLength             =  0;
  fDeltaE             = -7;
  fTotalStoppingPower = -8;
  fTotalNumberOfSteps = 0;
  fNumberOfSteps      = 0;
  fDepositedEnergy    = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{

 //*** WARNING: this line will kill all secondary e- and gammas ***

 if (    step->GetTrack()->GetDefinition()==G4Electron::ElectronDefinition()
      && step->GetTrack()->GetTrackID() != 1 ) 
    step->GetTrack()->SetTrackStatus(fStopAndKill);

 if (    step->GetTrack()->GetDefinition()==G4Gamma::GammaDefinition() ) 
    step->GetTrack()->SetTrackStatus(fStopAndKill);

 //
 
 G4double charge   = step->GetTrack()->GetDefinition()->GetPDGCharge();
 G4double steplen  = step->GetStepLength();
 G4double edep     =  step->GetTotalEnergyDeposit();
 G4int maxNumberOfSteps = -1;

 G4DNAGenericIonsManager *instance;
 instance = G4DNAGenericIonsManager::Instance();
   
 if ( 
 
      // ELECTRONS 
      (
        step->GetTrack()->GetTrackID() == 1 &&
        step->GetTrack()->GetParticleDefinition()==G4Electron::ElectronDefinition()
      )
     
      ||
     
      // PROTONS, HYDROGEN, ALPHA PARTICLES AND CHARGED STATES
      (
        charge != -1. &&
        step->GetTrack()->GetParticleDefinition()!=G4Gamma::GammaDefinition()       
      )
     
    )
 { 
     fSumOfStepLength = fSumOfStepLength + steplen;
     fDepositedEnergy = fDepositedEnergy + edep;
     fNumberOfSteps = fNumberOfSteps + 1;

     //G4cout << fSumOfStepLength << " " << fDepositedEnergy << " " << fNumberOfSteps << G4endl;
          
     // ELECTRONS
     
     if (step->GetTrack()->GetParticleDefinition()==G4Electron::ElectronDefinition())
      maxNumberOfSteps = 1000;

     // PROTONS, HYDROGEN
     
     if 
     (
       step->GetTrack()->GetParticleDefinition()==G4Proton::ProtonDefinition()
       ||
       step->GetTrack()->GetParticleDefinition()==instance->GetIon("hydrogen")
     )
       maxNumberOfSteps = 1000;
     
     // He0, He+, He2+
     
     if 
     (
       step->GetTrack()->GetParticleDefinition()==instance->GetIon("alpha++")
       ||
       step->GetTrack()->GetParticleDefinition()==instance->GetIon("alpha+")
       ||
       step->GetTrack()->GetParticleDefinition()==instance->GetIon("helium")
     )
       maxNumberOfSteps = 10000;
      

     // All
     
     if (fNumberOfSteps > maxNumberOfSteps)
     { 
        fLength = fSumOfStepLength;
        
        fTotalStoppingPower = fDepositedEnergy/fLength;
        //G4cout << fTotalStoppingPower/(MeV/cm) << G4endl;
 
        G4RunManager::GetRunManager()->AbortEvent();
        fSumOfStepLength=0.;
        fDepositedEnergy=0;
        fNumberOfSteps = 0;
 
     };
 };

}

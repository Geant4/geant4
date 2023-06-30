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
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StackingAction.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

#include "Run.hh"
#include "RunAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction() : G4UserStackingAction(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::~StackingAction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
  //    keep primary particle
  // count secondary particles
  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  auto energy_keV = (G4float)(aTrack->GetKineticEnergy() / keV);
  if (run->GetIsPIXE()) {
        const G4double threshold = 0.9;  //*keV;
    //          const G4double threshold = 0.0; //*keV;
    if (aTrack->GetParentID() > 0 && energy_keV > threshold
        && aTrack->GetDefinition() == G4Gamma::Gamma())
    {
      auto mx = (G4float)aTrack->GetMomentumDirection().x();
      auto my = (G4float)aTrack->GetMomentumDirection().y();
      auto mz = (G4float)aTrack->GetMomentumDirection().z();

      //      G4cout << std::setprecision(6) << energy_keV << G4endl;
      //      auto x = (G4float)(aTrack->GetPosition().x() / um);
      //      auto y = (G4float)(aTrack->GetPosition().y() / um);
      //      auto z = (G4float)(aTrack->GetPosition().z() / um);

      ParticleInfo gammaInfo(energy_keV, mx, my, mz);
      // ParticleInfo gammaInfo(energy_keV, mx, my, mz, x, y, z);
      run->FillGammaAtCreation(gammaInfo);
    }
  }

  // stack or delete secondaries
  G4ClassificationOfNewTrack status = fUrgent;
  return status;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

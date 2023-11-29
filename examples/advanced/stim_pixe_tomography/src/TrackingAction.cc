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
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include "G4Proton.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

#include "Run.hh"
#include "RunAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction() : G4UserTrackingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::~TrackingAction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  auto energy_keV = (G4float)(aTrack->GetKineticEnergy() / keV);
  if (run->GetIsPIXE()) {
    // gammas getting out of object
    // parentID > 0 => secondary
    const G4double threshold = 0.9;  //* keV;
    if (energy_keV > threshold && aTrack->GetDefinition() == G4Gamma::Gamma()
        && aTrack->GetParentID() > 0)
    {
      auto mx = (G4float)aTrack->GetMomentumDirection().x();
      auto my = (G4float)aTrack->GetMomentumDirection().y();
      auto mz = (G4float)aTrack->GetMomentumDirection().z();
      //
      //      auto x = (G4float)(aTrack->GetPosition().x() / um);
      //      auto y = (G4float)(aTrack->GetPosition().y() / um);
      //      auto z = (G4float)(aTrack->GetPosition().z() / um);

      //      ParticleInfo gammaInfo(energy_keV, mx, my, mz, x, y, z);
      ParticleInfo gammaInfo(energy_keV, mx, my, mz);
      run->FillGammaAtExit(gammaInfo);
    }
  }
  else {
    //
    // primary protons (getting out of object)
    // parentID > 0 => secondary; parentID = 0 => primary
    if (energy_keV > 0.0 && aTrack->GetDefinition() == G4Proton::Proton()
        && aTrack->GetParentID() == 0) {
      auto mx = (G4float)aTrack->GetMomentumDirection().x();
      auto my = (G4float)aTrack->GetMomentumDirection().y();
      auto mz = (G4float)aTrack->GetMomentumDirection().z();

      //      auto x = (G4float)(aTrack->GetPosition().x() / um);
      //      auto y = (G4float)(aTrack->GetPosition().y() / um);
      //      auto z = (G4float)(aTrack->GetPosition().z() / um);

      //      G4cout << std::setprecision(6) << x << " " << y << " " << z
      //             << G4endl;

      //      ParticleInfo protonInfo(energy_keV, mx, my, mz, x, y, z);
      ParticleInfo protonInfo(energy_keV, mx, my, mz);
      run->FillProtonAtExit(protonInfo);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

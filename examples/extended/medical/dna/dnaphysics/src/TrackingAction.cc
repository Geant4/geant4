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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publications:
// Med. Phys. 45 (2018) e722-e739
// Phys. Med. 31 (2015) 861-874
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file TrackingAction.cc
/// \brief Implementation of the TrackingAction class

#include "TrackingAction.hh"

#include "G4Alpha.hh"
#include "G4AnalysisManager.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4PhysicalConstants.hh"
#include "G4Proton.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  G4double flagParticle = -1.;
  G4double x, y, z, dirx, diry, dirz;

  G4ParticleDefinition* partDef = aTrack->GetDynamicParticle()->GetDefinition();

  if (partDef == G4Gamma::GammaDefinition()) flagParticle = 0;

  if (partDef == G4Electron::ElectronDefinition()) flagParticle = 1;

  if (partDef == G4Proton::ProtonDefinition()) flagParticle = 2;

  if (partDef == G4Alpha::AlphaDefinition()) flagParticle = 4;

  G4DNAGenericIonsManager* instance;
  instance = G4DNAGenericIonsManager::Instance();

  if (partDef == instance->GetIon("hydrogen")) flagParticle = 3;

  if (partDef == instance->GetIon("alpha+")) flagParticle = 5;

  if (partDef == instance->GetIon("helium")) flagParticle = 6;

  //

  x = aTrack->GetPosition().x() / nanometer;
  y = aTrack->GetPosition().y() / nanometer;
  z = aTrack->GetPosition().z() / nanometer;

  dirx = aTrack->GetMomentumDirection().x();
  diry = aTrack->GetMomentumDirection().y();
  dirz = aTrack->GetMomentumDirection().z();

  // Call analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Fill track information ntuple
  analysisManager->FillNtupleDColumn(1, 0, flagParticle);
  analysisManager->FillNtupleDColumn(1, 1, x);
  analysisManager->FillNtupleDColumn(1, 2, y);
  analysisManager->FillNtupleDColumn(1, 3, z);
  analysisManager->FillNtupleDColumn(1, 4, dirx);
  analysisManager->FillNtupleDColumn(1, 5, diry);
  analysisManager->FillNtupleDColumn(1, 6, dirz);
  analysisManager->FillNtupleDColumn(1, 7, aTrack->GetKineticEnergy() / eV);
  analysisManager->FillNtupleIColumn(1, 8, aTrack->GetTrackID());
  analysisManager->FillNtupleIColumn(1, 9, aTrack->GetParentID());
  analysisManager->AddNtupleRow(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track*) {}

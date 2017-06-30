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
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// and papers
// M. Batmunkh et al. J Radiat Res Appl Sci 8 (2015) 498-507
// O. Belov et al. Physica Medica 32 (2016) 1510-1520
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// -------------------------------------------------------------------
// November 2016
// -------------------------------------------------------------------
//
// $ID$
/// \file TrackingAction.cc 
/// \brief Implementation of the TrackingAction class

#include "TrackingAction.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "G4Region.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4RegionStore.hh"
//
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"  
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
//
// Bosons
#include "G4ChargedGeantino.hh"
#include "G4Geantino.hh"
#include "G4Gamma.hh"
// leptons
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"
// Mesons
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"
#include "G4Eta.hh"
#include "G4EtaPrime.hh"

#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonZero.hh"
#include "G4AntiKaonZero.hh"
#include "G4KaonZeroLong.hh"
#include "G4KaonZeroShort.hh"
// Baryons
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Neutron.hh"
#include "G4AntiNeutron.hh"
// Nuclei
#include "G4Alpha.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4GenericIon.hh"
#include "G4Step.hh"
#include "G4StepStatus.hh"
#include "Run.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction() :
    G4UserTrackingAction(),
    RunInitObserver(), fpTargetRegion(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::~TrackingAction()
{
    fpTargetRegion = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
 
// Target volumes
  G4VPhysicalVolume* volumeStep = track->GetTouchableHandle()->GetVolume();
  G4VPhysicalVolume* volumeMedium = 
           G4PhysicalVolumeStore::GetInstance()->GetVolume("Medium");
  G4VPhysicalVolume* volumeSlice = 
           G4PhysicalVolumeStore::GetInstance()->GetVolume("BoundingSlice"); 
  Run* run = static_cast<Run*>(
           G4RunManager::GetRunManager()->GetNonConstCurrentRun());   
 
    //count secondary particles
    if (track->GetTrackID() == 1) return;  
    G4String name   = track->GetDefinition()->GetParticleName();
    G4double energy = track->GetKineticEnergy();
 
  // outside bounding slice
     //if (volumeStep == volumeMedium) run->ParticleCount(name,energy);
  // inside bounding slice 
     //if (volumeStep != volumeMedium) run->ParticleCountNeuron(name,energy);
    
 // particles outside neuron structure
 if (volumeStep == volumeMedium || volumeStep == volumeSlice)
 {
  run->ParticleCount(name,energy);
 }
 // count secondary particles in neuron
 else //if (volumeStep != volumeMedium && volumeStep != volumeSlice)
 {
  run->ParticleCountNeuron(name,energy);
 }

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{

    if(track->GetTrackID() == 1)
    {
        G4cout<<"End of tracking primary particle, its final energy is :"
            << G4BestUnit(track->GetKineticEnergy(), "Energy")<< G4endl;
    }
 
  Run* run = static_cast<Run*>(
  G4RunManager::GetRunManager()->GetNonConstCurrentRun());   
  G4double tracklen = track->GetTrackLength();
  if (track->GetTrackID() == 1) run->SetTrackLength(tracklen); 
  //G4cout<< " track length = " << tracklen / um << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void
TrackingAction::Initialize()
{
  //fpTargetRegion = G4RegionStore::GetInstance()->GetRegion("BoundingSlice");
}

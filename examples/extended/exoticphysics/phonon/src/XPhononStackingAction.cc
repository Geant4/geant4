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
/// \file exoticphysics/phonon/src/XPhononStackingAction.cc
/// \brief Implementation of the XPhononStackingAction class
///     This stacking action is necessary to ensure that velocity and 
///     propagation direction are set properly for phonons created with
///     G4ParticleGun
//
//

#include "XPhononStackingAction.hh"
#include "G4LatticeManager.hh"
#include "G4PhononLong.hh"
#include "G4PhononPolarization.hh"
#include "G4PhononTrackMap.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XPhononStackingAction::XPhononStackingAction() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XPhononStackingAction::~XPhononStackingAction() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ClassificationOfNewTrack 
XPhononStackingAction::ClassifyNewTrack(const G4Track* aTrack) {
  G4ClassificationOfNewTrack classification = fUrgent;

  if (aTrack->GetParentID() == 0) {
    //Obtain LatticeManager for phonon dynamics
    G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
    
    G4int pol = G4PhononPolarization::Get(aTrack->GetDefinition());
    
    //Compute random wave-vector (override whatever ParticleGun did)
    G4ThreeVector Ran = G4RandomDirection();
    
    //Store wave-vector as track information
    G4PhononTrackMap* theKmap = G4PhononTrackMap::GetPhononTrackMap();
    theKmap->SetK(aTrack, Ran);
    
    //Compute direction of propagation from wave vector
    G4ThreeVector momentumDir = LM->MapKtoVDir(aTrack->GetVolume(), pol, Ran);
    
    //Compute true velocity of propagation
    G4double velocity = LM->MapKtoV(aTrack->GetVolume(), pol, Ran);
    
    //cast to non-const pointer so we can set the velocity
    G4Track* theTrack = const_cast<G4Track*>(aTrack);

    theTrack->SetMomentumDirection(momentumDir);
    theTrack->SetVelocity(velocity);
    theTrack->UseGivenVelocity(true);
  }
  
  return classification; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


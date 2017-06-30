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
/// \file processes/phonon/src/G4VPhononProcess.cc
/// \brief Implementation of the G4VPhononProcess base class
//
// $Id: G4VPhononProcess.cc 76777 2013-11-15 16:20:56Z mkelsey $
//
// 20131111  Add verbosity to report creating secondaries

#include "G4VPhononProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4ExceptionSeverity.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhononLong.hh"
#include "G4PhononPolarization.hh"
#include "G4PhononTrackMap.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4ProcessType.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"

namespace {
  const G4ThreeVector nullVec(0.,0.,0.);	// For convenience below
}

// Constructor and destructor

G4VPhononProcess::G4VPhononProcess(const G4String& processName)
  : G4VDiscreteProcess(processName, fPhonon),
    trackKmap(G4PhononTrackMap::GetInstance()), theLattice(0),
    currentTrack(0) {}

G4VPhononProcess::~G4VPhononProcess() {;}


// Only applies to the known phonon polarization states

G4bool G4VPhononProcess::IsApplicable(const G4ParticleDefinition& aPD) {
  return (&aPD==G4PhononLong::Definition() ||
	  &aPD==G4PhononTransFast::Definition() ||
	  &aPD==G4PhononTransSlow::Definition() );
}


// Initialize wave vectors for currently active track(s)

void G4VPhononProcess::StartTracking(G4Track* track) {
  G4VProcess::StartTracking(track);	// Apply base class actions

  // FIXME:  THE WAVEVECTOR SHOULD BE COMPUTED BY INVERTING THE K/V MAP
  if (!trackKmap->Find(track)) 
    trackKmap->SetK(track, track->GetMomentumDirection());

  currentTrack = track;			// Save for use by EndTracking

  // Fetch lattice for current track once, use in subsequent steps
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  theLattice = LM->GetLattice(track->GetVolume());
}

void G4VPhononProcess::EndTracking() {
  G4VProcess::EndTracking();		// Apply base class actions
  trackKmap->RemoveTrack(currentTrack);
  currentTrack = 0;
  theLattice = 0;
}


// For convenience, map phonon type to polarization code

G4int G4VPhononProcess::GetPolarization(const G4Track& track) const {
  return G4PhononPolarization::Get(track.GetParticleDefinition());
}


// Generate random polarization from density of states

G4int G4VPhononProcess::ChoosePolarization(G4double Ldos, G4double STdos,
					   G4double FTdos) const {
  G4double norm = Ldos + STdos + FTdos;
  G4double cProbST = STdos/norm;
  G4double cProbFT = FTdos/norm + cProbST;

  // NOTE:  Order of selection done to match previous random sequences
  G4double modeMixer = G4UniformRand();
  if (modeMixer<cProbST) return G4PhononPolarization::TransSlow;
  if (modeMixer<cProbFT) return G4PhononPolarization::TransFast;
  return G4PhononPolarization::Long;
}


// Create new secondary track from phonon configuration

G4Track* G4VPhononProcess::CreateSecondary(G4int polarization,
					   const G4ThreeVector& waveVec,
					   G4double energy) const {
  if (verboseLevel>1) {
    G4cout << GetProcessName() << " CreateSecondary pol " << polarization
	   << " K " << waveVec << " E " << energy << G4endl;
  }

  G4ThreeVector vgroup = theLattice->MapKtoVDir(polarization, waveVec);
  if (verboseLevel>1) G4cout << " MapKtoVDir returned " << vgroup << G4endl;

  vgroup = theLattice->RotateToGlobal(vgroup);
  if (verboseLevel>1) G4cout << " RotateToGlobal returned " << vgroup << G4endl;

  if (verboseLevel && std::fabs(vgroup.mag()-1.) > 0.01) {
    G4cout << "WARNING: " << GetProcessName() << " vgroup not a unit vector: "
	   << vgroup << G4endl;
  }

  G4ParticleDefinition* thePhonon = G4PhononPolarization::Get(polarization);

  // Secondaries are created at the current track coordinates
  G4Track* sec = new G4Track(new G4DynamicParticle(thePhonon, vgroup, energy),
			     currentTrack->GetGlobalTime(),
			     currentTrack->GetPosition());

  // Store wavevector in lookup table for future tracking
  trackKmap->SetK(sec, theLattice->RotateToGlobal(waveVec));

  if (verboseLevel>1) {
    G4cout << GetProcessName() << " secondary K rotated to "
	   << trackKmap->GetK(sec) << G4endl;
  }

  sec->SetVelocity(theLattice->MapKtoV(polarization, waveVec));    
  sec->UseGivenVelocity(true);

  return sec;
}

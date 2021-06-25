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
// --------------------------------------------------------------------
// GEANT4 class header file
//
// Class Description:
// Base class for a clipboard for communicating between quantum entangled tracks.
//
// ------------------ G4VEntanglementClipBoard ------------------
//
// Author: J.Allison, May 2017
//
// --------------------------------------------------------------------

// Usage:
//
// In the method that generates entangled secondaries
// (See, for example, G4eplusAnnihilation::AtRestDoIt)
// Make a shared pointer to a clip board and attach it to the tracks through
// G4EntanglementAuxInfo. That way the clip board lasts the life of both tracks.
// (G4XXXEntanglementClipBoard is your class inherited from this class)
// (See, for example, G4eplusAnnihilationEntanglementClipBoard)
//  auto clipBoard = std::make_shared<G4XXXEntanglementClipBoard>();
// For each secondary
//  G4Track* track = new G4Track(...
//  ...
//  clipBoard->SetTrackA(track);
//  track->SetAuxiliaryTrackInformation(0,new G4EntanglementAuxInfo(clipBoard));
// Then repeat for track B
//
// In the method that does the quantum "measurement"
// (See, for example, G4LivermorePolarizedComptonModel::SampleSecondaries)
//  const auto* auxInfo = fParticleChange->GetCurrentTrack()->GetAuxiliaryTrackInformation(0);
//  if (auxInfo) {
//    const auto* entanglementAuxInfo = dynamic_cast<const G4EntanglementAuxInfo*>(auxInfo);
//    if (entanglementAuxInfo) {
//      auto* clipBoard = dynamic_cast<G4XXXEntanglementClipBoard*>
//      (entanglementAuxInfo->GetEntanglementClipBoard());
//      if (clipBoard) {
//        if (clipBoard->IsTrack1Measurement()) {
//          if (clipBoard->GetTrackB() == fParticleChange->GetCurrentTrack()) {
//            clipBoard->ResetTrack1Measurement();
//            // Store information
//            ...
//          }
//        } else if (clipBoard->IsTrack2Measurement()) {
//          if (clipBoard->GetTrackA() == fParticleChange->GetCurrentTrack()) {
//            clipBoard->ResetTrack2Measurement();
//            // Retrieve information and apply quantum mechanics
//            ...
//          }
//        }
//      }
//    }
//  }

#ifndef G4VEntanglementClipBoard_hh
#define G4VEntanglementClipBoard_hh

#include "globals.hh"

class G4ParticleDefinition;
class G4Track;

class G4VEntanglementClipBoard {

public:
  
  G4VEntanglementClipBoard()
  : fpParentParticleDefinition(0)
  , fTrackA(0)
  , fTrackB(0)
  , fTrack1Measurement(true)
  , fTrack2Measurement(true)
  {}
  virtual ~G4VEntanglementClipBoard() {}

  void SetParentParticleDefinition(G4ParticleDefinition* p)
  {fpParentParticleDefinition = p;}
  G4ParticleDefinition* GetParentParticleDefinition() const
  {return fpParentParticleDefinition;}

  void SetTrackA(const G4Track* track) {fTrackA = track;}
  void SetTrackB(const G4Track* track) {fTrackB = track;}
  const G4Track* GetTrackA() const {return fTrackA;}
  const G4Track* GetTrackB() const {return fTrackB;}

  // The entanglement-sensitive process is responsible for setting this.
  void ResetTrack1Measurement() {fTrack1Measurement = false;}
  void ResetTrack2Measurement() {fTrack2Measurement = false;}
  G4bool IsTrack1Measurement() const {return fTrack1Measurement;}
  G4bool IsTrack2Measurement() const {return fTrack2Measurement;}
  
private:

  G4ParticleDefinition* fpParentParticleDefinition;

  // Use pointer values to identify tracks. We don't know in which order the
  // tracks will be processed so let us call them A & B.
  const G4Track* fTrackA;
  const G4Track* fTrackB;

  // False until first measurement...
  G4bool fTrack1Measurement;  // ...of the first encountered track
  G4bool fTrack2Measurement;  // ...of the second encountered track
  
};

#endif

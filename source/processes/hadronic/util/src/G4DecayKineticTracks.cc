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
// Utility class to process contents of G4KineticTrackVector (input to
// models' ::Propagate() interface) and decay any short-lived resonances.
// Resulting daughters are added to vector (no function or return needed).
//
// Author:  Michael Kelsey <kelsey@slac.stanford.edu>

#include "G4DecayKineticTracks.hh"
#include "G4KineticTrackVector.hh"
#include "G4KineticTrack.hh"


// Decay all input tracks, put daughters onto end of list

G4DecayKineticTracks::G4DecayKineticTracks(G4KineticTrackVector *tracks) {

  if (tracks) Decay(tracks);
}

void G4DecayKineticTracks::Decay(G4KineticTrackVector *tracks) const {

  if (!tracks) return;

  G4KineticTrackVector* daughters = 0;
  for (size_t i=0; i<tracks->size(); ++i) {
    G4KineticTrack* track = (*tracks)[i];
    if (!track) continue;

    // Select decay of current track, put daughters at end of vector
    daughters = track->GetDefinition()->IsShortLived() ? track->Decay() : 0;
  
    if (daughters) {
      // Use the integer round mass in keV to get an unique ID for the parent resonance
      G4int uniqueID = static_cast< G4int >( round( track->Get4Momentum().mag() / CLHEP::keV ) );
 
      // Assign to the daughters the creator model ID of their parent
      for (size_t k=0; k<daughters->size(); ++k) {
	G4KineticTrack* aDaughter = (*daughters)[k];
	if (aDaughter) {
          aDaughter->SetCreatorModelID(track->GetCreatorModelID());
          aDaughter->SetParentResonanceDef(track->GetDefinition());
          aDaughter->SetParentResonanceID(uniqueID);          
        }
      }
      
      tracks->insert(tracks->end(), daughters->begin(), daughters->end());
      delete track;		// Remove parent track
      delete daughters;
      (*tracks)[i] = nullptr;	// Flag parent's slot for removal
    }
  }

  // Find and remove null pointers created by decays above
  for (G4int j=(G4int)tracks->size()-1; j>=0; --j) {
    if (nullptr == (*tracks)[j]) tracks->erase(tracks->begin()+j);
  }
}

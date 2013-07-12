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

#include "G4DecKineticTracks.hh"
#include "G4KineticTrackVector.hh"
#include "G4KineticTrack.hh"


// Decay all input tracks, put daughters onto end of list

G4DecKineticTracks::G4DecKineticTracks(G4KineticTrackVector *tracks) {
  if (tracks) Decay(tracks);
}

void G4DecKineticTracks::Decay(G4KineticTrackVector *tracks) const {
  if (!tracks) return;

/*
G4cout<<"Initial particles"<<G4endl;
for (size_t i=0; i<tracks->size(); ++i) 
{
 G4KineticTrack* track = (*tracks)[i];
 G4cout<<i<<" "<<track->GetDefinition()->GetParticleName()<<" "<<track->GetDefinition()->IsShortLived()<<G4endl;
}
*/

  G4KineticTrackVector* daughters = 0;
  for (size_t i=0; i<tracks->size(); ++i) {
    G4KineticTrack* track = (*tracks)[i];
    if (!track) continue;

    // Select decay of current track, put daughters at end of vector
/*
G4cout<<" ----- "
<<track->GetDefinition()->GetParticleName()<<" "
<<track->GetDefinition()->IsShortLived()   <<" "
<<track->GetDefinition()->GetPDGMass()     <<" "
<<track->GetDefinition()->GetPDGCharge()   <<G4endl;
*/
    daughters = track->GetDefinition()->IsShortLived() ? track->Decay() : 0;
if(track->GetDefinition()->GetParticleName() == "eta"      ) daughters = track->Decay();
if(track->GetDefinition()->GetParticleName() == "eta_prime") daughters = track->Decay();
    if (daughters) {
/*
G4cout<<"Dauthers"<<G4endl;
for (size_t j=0; j<daughters->size(); ++j) 
{

 G4cout<<j<<" *** "<<(*daughters)[j]->GetDefinition()->GetParticleName()<<G4endl;
}
*/
      tracks->insert(tracks->end(), daughters->begin(), daughters->end());
      delete track;		// Remove parent track
      delete daughters;
      (*tracks)[i] = NULL;	// Flag parent's slot for removal
    }
  }

  // Find and remove null pointers created by decays above
  for (int j=tracks->size()-1; j>=0; --j) {
    if (NULL == (*tracks)[j]) tracks->erase(tracks->begin()+j);
  }
}

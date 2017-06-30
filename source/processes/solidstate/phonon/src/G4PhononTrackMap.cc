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
/// \file processes/phonon/src/G4PhononTrackMap.hh
/// \brief Implementation of the G4PhononTrackMap base class
//
// $Id: G4PhononTrackMap.cc 76492 2013-11-11 17:15:04Z mkelsey $
//
// 20131111  Move Clear() function to .cc file

#include "G4PhononTrackMap.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include <map>

namespace { const G4ThreeVector nullVec(0.,0.,0.); }	// For convenience below

// Singleton instance, must be created at first request

G4ThreadLocal G4PhononTrackMap* G4PhononTrackMap::theTrackMap = 0;


// Pseudo-constructor creates singleton instance

G4PhononTrackMap* G4PhononTrackMap::GetPhononTrackMap() {
  if (!theTrackMap) theTrackMap = new G4PhononTrackMap;
  return theTrackMap;
}

void G4PhononTrackMap::Clear() {
  theMap.clear();			// Remove all entries from map
}



// Check if specified track is already loaded

G4bool G4PhononTrackMap::Find(const G4Track* track) const {
  return (!track || theMap.find(track) != theMap.end());
}


// Remove specified track from map (used by EndTracking)

void G4PhononTrackMap::RemoveTrack(const G4Track* track) {
  TrkIDKmap::iterator entry = theMap.find(track);
  if (entry != theMap.end()) theMap.erase(entry);
}


// Update the wavevector for specified track, add track if non-existent

void G4PhononTrackMap::SetK(const G4Track* track, const G4ThreeVector& K) {
  if (track) theMap[track] = K;
}


// Access current wavevector for specified track (NULL if doesn't exist)

const G4ThreeVector& G4PhononTrackMap::GetK(const G4Track* track) const {
  TrkIDKmap::const_iterator entry = theMap.find(track);
  return (entry != theMap.end()) ? entry->second : nullVec;
}
  

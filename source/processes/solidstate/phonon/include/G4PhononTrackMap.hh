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
/// \file processes/phonon/include/G4PhononTrackMap.hh
/// \brief Definition of the G4PhononTrackMap base class
//
// $Id: G4PhononTrackMap.hh 76492 2013-11-11 17:15:04Z mkelsey $
//
// 20131111  Move implementation of Clear() to .cc file

#ifndef G4PhononTrackMap_h
#define G4PhononTrackMap_h 1

#include "G4ThreeVector.hh"
#include <map>

class G4Track;

class G4PhononTrackMap {
public:
  typedef std::map<const G4Track*, G4ThreeVector> TrkIDKmap;
  static G4ThreadLocal G4PhononTrackMap* theTrackMap;

public:
  static G4PhononTrackMap* GetPhononTrackMap();	// Synonyms for access
  static G4PhononTrackMap* GetInstance() { return GetPhononTrackMap(); }

  // Update the wavevector for specified track, add track if not found
  void SetK(const G4Track* track, const G4ThreeVector& K);
  void SetK(const G4Track& track, const G4ThreeVector& K) { SetK(&track, K); }

  // Access current wavevector for specified track (NULL if doesn't exist)
  const G4ThreeVector& GetK(const G4Track* track) const;
  const G4ThreeVector& GetK(const G4Track& track) const { return GetK(&track); }

  // Check if specified track is already loaded
  G4bool Find(const G4Track* track) const;
  G4bool Find(const G4Track& track) const { return Find(&track); }

  // Remove specified track from map (used by EndTracking)
  void RemoveTrack(const G4Track* track);

  void Clear();			// Remove all entries from map

private:
  TrkIDKmap theMap;		// Associate track ID numbers with vectors

private:
  G4PhononTrackMap() { Clear(); }		// Ensure map is empty
  ~G4PhononTrackMap() {;}
};
#endif	/* G4PhononTrackMap_h */

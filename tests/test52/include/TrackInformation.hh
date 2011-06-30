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
#ifndef TRACKINFORMATION_HH
#define TRACKINFORMATION_HH

#include "G4VUserTrackInformation.hh"
#include "G4Allocator.hh"
#include "globals.hh"


class TrackInformation : public G4VUserTrackInformation {
 
 public:
   TrackInformation() : hitTarget(false), escapeEnergy(0) {}
   ~TrackInformation() {}
   inline void *operator new(size_t);
   inline void operator delete(void* trackInfo);

   void Print() const {}

 private:
   G4bool hitTarget;
   G4double escapeEnergy;

 public:
   inline void SetHitTarget() { 
     hitTarget = true;
   }
   inline G4bool HitTarget() const {
     return hitTarget;
   }
   inline void SetEscapeEnergy(G4double en) {
     escapeEnergy = en;
   }
   inline G4double EscapeEnergy() {
     return escapeEnergy;
   }
};

extern G4Allocator<TrackInformation> TrackInformationAllocator;

inline void* TrackInformation::operator new(size_t) { 
void* trackInfo;
  trackInfo = (void*) TrackInformationAllocator.MallocSingle();
  return trackInfo;
}

inline void TrackInformation::operator delete(void* trackInfo) { 
  TrackInformationAllocator.FreeSingle((TrackInformation*) trackInfo);
}

#endif // TRACKINFORMATION_HH

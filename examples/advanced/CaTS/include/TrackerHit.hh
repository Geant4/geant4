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

// ********************************************************************
//
//  CaTS (Calorimetry and Tracking Simulation)
//
//  Authors : Hans Wenzel
//            Soon Yung Jun
//            (Fermi National Accelerator Laboratory)
//
// History
//   October 18th, 2021 : first implementation
//
// ********************************************************************
//
/// \file TrackerHit.hh
/// \brief Definition of the CaTS::TrackerHit class

#pragma once

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class TrackerHit : public G4VHit
{
 public:
  TrackerHit();
  ~TrackerHit() = default;
  TrackerHit(const TrackerHit&);
  const TrackerHit& operator=(const TrackerHit&);
  G4bool operator==(const TrackerHit&) const;
  inline void* operator new(size_t);
  inline void operator delete(void*);
  void Draw() final;

  inline void Print() final
  {
    G4cout << "TrackerHit  id:  " << fid << " Edep: " << fEdep
           << " X: " << fposition.getX() << " Y: " << fposition.getY()
           << " Z: " << fposition.getZ() << " time: " << ftime << G4endl;
  }

  TrackerHit(G4double edep, G4ThreeVector position, G4double time);
  inline void SetEdep(G4double Edep) { fEdep = Edep; }
  inline G4double GetEdep() { return fEdep; }
  inline void SetTime(G4double time) { ftime = time; }
  inline G4double GetTime() const { return ftime; }
  inline void SetPosition(G4ThreeVector position) { fposition = position; }
  inline G4ThreeVector GetPosition() const { return fposition; }

 private:
  G4int fid{ 0 };
  G4double fEdep{ 0 };
  G4ThreeVector fposition{ 0, 0, 0 };
  G4double ftime{ 0 };
};

using TrackerHitsCollection = G4THitsCollection<TrackerHit>;
extern G4ThreadLocal G4Allocator<TrackerHit>* TrackerHitAllocator;

inline void* TrackerHit::operator new(size_t)
{
  if(!TrackerHitAllocator)
  {
    TrackerHitAllocator = new G4Allocator<TrackerHit>;
  }
  return (void*) TrackerHitAllocator->MallocSingle();
}

inline void TrackerHit::operator delete(void* aHit)
{
  TrackerHitAllocator->FreeSingle((TrackerHit*) aHit);
}

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
#ifndef TARGETLAYERHIT_HH
#define TARGETLAYERHIT_HH

#include <vector>
#include <utility>
#include "G4VHit.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

typedef std::pair<G4String,G4double> NewParticle;
typedef std::vector<NewParticle> Secondaries;


class TargetLayerHit : public G4VHit {

 public:
   TargetLayerHit();
   TargetLayerHit(const TargetLayerHit& right);
   ~TargetLayerHit();
   
   const TargetLayerHit& operator= (const TargetLayerHit& right);
   int operator==(const TargetLayerHit& right) const;

   inline void* operator new(size_t);
   inline void operator delete(void* hit);

 private:
   G4String volName;
   G4double energyDeposit;
   G4ThreeVector positionPre;
   G4ThreeVector positionPost;
   G4String particle;
   G4int trackID;
   Secondaries secParticles;

 public:

   void Draw();
   void Print() {}

   inline void SetVolumeName(G4String name) {
     volName = name;
   }
   inline void StoreEnergyDeposit(G4double energy) {
     energyDeposit = energy;
   }
   inline void StorePositionPre(G4ThreeVector pos) {
     positionPre = pos;
   }
   inline void StorePositionPost(G4ThreeVector pos) {
     positionPost = pos;
   }
   inline void StoreParticleType(G4String par) {
     particle = par;
   }
   inline void StoreTrackID(G4int id) {
     trackID = id;
   }
   inline void StoreSecParticle(G4String name, G4double energy) {
     secParticles.push_back(std::make_pair(name,energy)); 
   }
   inline G4double GetEnergyDeposit() {
     return energyDeposit;
   }
   inline G4ThreeVector GetPosition() {
     return positionPre;
   }
   inline G4double GetXCoordPre() {
     return positionPre.x();
   }
   inline G4double GetYCoordPre() {
     return positionPre.y();
   }
   inline G4double GetZCoordPre() {
     return positionPre.z();
   }
   inline G4double GetXCoordPost() {
     return positionPost.x();
   }
   inline G4double GetYCoordPost() {
     return positionPost.y();
   }
   inline G4double GetZCoordPost() {
     return positionPost.z();
   }
   inline G4String GetParticleType() {
     return particle;
   }
   inline G4int GetTrackID() {
     return trackID; 
   }
   inline Secondaries& GetSecParticles() {
     return secParticles;
   }
};


extern G4Allocator<TargetLayerHit> TargetLayerHitAllocator;

inline void* TargetLayerHit::operator new(size_t) {

  void* hit;
  hit = (void*) TargetLayerHitAllocator.MallocSingle();
  return hit;
}

inline void TargetLayerHit::operator delete(void* hit) {

  TargetLayerHitAllocator.FreeSingle((TargetLayerHit*) hit);
}

#endif // TARGETLAYERHIT_HH

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
/// \file optical/LXe/include/LXePMTHit.hh
/// \brief Definition of the LXePMTHit class
//
//
#ifndef LXePMTHit_h
#define LXePMTHit_h 1

#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"
#include "G4THitsCollection.hh"
#include "G4VHit.hh"
#include "G4VPhysicalVolume.hh"

class LXePMTHit : public G4VHit
{
 public:
  LXePMTHit() = default;
  LXePMTHit(const LXePMTHit& right);
  ~LXePMTHit() override = default;

  const LXePMTHit& operator=(const LXePMTHit& right);
  G4bool operator==(const LXePMTHit& right) const;

  inline void* operator new(size_t);
  inline void operator delete(void* aHit);

  void Draw() override;
  void Print() override;

  inline void SetDrawit(G4bool b) { fDrawit = b; }
  inline G4bool GetDrawit() { return fDrawit; }

  inline void IncPhotonCount() { ++fPhotons; }
  inline G4int GetPhotonCount() { return fPhotons; }

  inline void SetPMTNumber(G4int n) { fPmtNumber = n; }
  inline G4int GetPMTNumber() { return fPmtNumber; }

  inline void SetPMTPhysVol(G4VPhysicalVolume* physVol)
  {
    fPhysVol = physVol;
  }
  inline G4VPhysicalVolume* GetPMTPhysVol() { return fPhysVol; }

  inline void SetPMTPos(G4double x, G4double y, G4double z)
  {
    fPos = G4ThreeVector(x, y, z);
  }

  inline G4ThreeVector GetPMTPos() { return fPos; }

 private:
  G4int fPmtNumber = -1;
  G4int fPhotons = 0;
  G4ThreeVector fPos;
  G4VPhysicalVolume* fPhysVol = nullptr;
  G4bool fDrawit = false;
};

typedef G4THitsCollection<LXePMTHit> LXePMTHitsCollection;

extern G4ThreadLocal G4Allocator<LXePMTHit>* LXePMTHitAllocator;

inline void* LXePMTHit::operator new(size_t)
{
  if(!LXePMTHitAllocator)
    LXePMTHitAllocator = new G4Allocator<LXePMTHit>;
  return (void*) LXePMTHitAllocator->MallocSingle();
}

inline void LXePMTHit::operator delete(void* aHit)
{
  LXePMTHitAllocator->FreeSingle((LXePMTHit*) aHit);
}

#endif

//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//

#ifndef LXeScintHit_h
#define LXeScintHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"

class LXeScintHit : public G4VHit
{
public:
  
  LXeScintHit();
  LXeScintHit(G4VPhysicalVolume* pVol);
  ~LXeScintHit();
  LXeScintHit(const LXeScintHit &right);
  const LXeScintHit& operator=(const LXeScintHit &right);
  G4int operator==(const LXeScintHit &right) const;

  inline void *operator new(size_t);
  inline void operator delete(void *aHit);
  
  void Draw();
  void Print();

  inline void SetEdep(G4double de) { edep = de; }  
  inline void AddEdep(G4double de) { edep += de; }
  inline G4double GetEdep() { return edep; }
  
  inline void SetPos(G4ThreeVector xyz) { pos = xyz; }
  inline G4ThreeVector GetPos() { return pos; }

  inline const G4VPhysicalVolume * GetPhysV() { return physVol; }

private:
  G4double edep;
  G4ThreeVector pos;
  const G4VPhysicalVolume* physVol;

};

typedef G4THitsCollection<LXeScintHit> LXeScintHitsCollection;

extern G4Allocator<LXeScintHit> LXeScintHitAllocator;

inline void* LXeScintHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) LXeScintHitAllocator.MallocSingle();
  return aHit;
}

inline void LXeScintHit::operator delete(void *aHit)
{
  LXeScintHitAllocator.FreeSingle((LXeScintHit*) aHit);
}

#endif




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




#ifndef ExN04MuonHit_h
#define ExN04MuonHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class ExN04MuonHit : public G4VHit
{
  public:

      ExN04MuonHit();
      ~ExN04MuonHit();
      ExN04MuonHit(const ExN04MuonHit &right);
      const ExN04MuonHit& operator=(const ExN04MuonHit &right);
      int operator==(const ExN04MuonHit &right) const;


      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      void Draw();
      void Print();

  private:
      G4double edep;
      G4ThreeVector pos;

  public:
      inline void SetEdep(G4double de)
      { edep = de; }
      inline void AddEdep(G4double de)
      { edep += de; }
      inline G4double GetEdep()
      { return edep; }
      inline void SetPos(G4ThreeVector xyz)
      { pos = xyz; }
      inline G4ThreeVector GetPos()
      { return pos; }

};

typedef G4THitsCollection<ExN04MuonHit> ExN04MuonHitsCollection;

extern G4Allocator<ExN04MuonHit> ExN04MuonHitAllocator;

inline void* ExN04MuonHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) ExN04MuonHitAllocator.MallocSingle();
  return aHit;
}

inline void ExN04MuonHit::operator delete(void *aHit)
{
  ExN04MuonHitAllocator.FreeSingle((ExN04MuonHit*) aHit);
}

#endif



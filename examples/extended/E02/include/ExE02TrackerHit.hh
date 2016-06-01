
#ifndef ExE02TrackerHit_h
#define ExE02TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh" 
#include "G4ThreeVector.hh"

class ExE02TrackerHit : public G4VHit
{
  public:

      ExE02TrackerHit();
      ~ExE02TrackerHit();
      ExE02TrackerHit(const ExE02TrackerHit &right);
      const ExE02TrackerHit& operator=(const ExE02TrackerHit &right);
      int operator==(const ExE02TrackerHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      void Draw();
      void Print();

  private:
      G4double edep;
      G4ThreeVector pos;

  public:
      inline void SetEdep(G4double de)
      { edep = de; };
      inline G4double GetEdep()
      { return edep; };
      inline void SetPos(G4ThreeVector xyz)
      { pos = xyz; };
      inline G4ThreeVector GetPos()
      { return pos; };

};

typedef G4THitsCollection<ExE02TrackerHit> ExE02TrackerHitsCollection;

extern G4Allocator<ExE02TrackerHit> ExE02TrackerHitAllocator;

inline void* ExE02TrackerHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) ExE02TrackerHitAllocator.MallocSingle();
  return aHit;
}

inline void ExE02TrackerHit::operator delete(void *aHit)
{
  ExE02TrackerHitAllocator.FreeSingle((ExE02TrackerHit*) aHit);
}

#endif



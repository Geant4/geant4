#ifndef MLHit_h
#define MLHit_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include <vector>

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

#include "G4HCofThisEvent.hh"

typedef std::vector<G4ThreeVector> Fluxhit;
////////////////////////////////////////////////////////////////////////////////
//
class MLHit : public G4VHit
{
public:
  MLHit ();
  ~MLHit ();
  MLHit (const MLHit&);
  const MLHit& operator= (const MLHit&);
  int operator== (const MLHit&) const;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  void SetEdep(G4int i, G4double de, G4double w) {Edep = de; Weight = w; Slayer = i;};
  G4double GetEdep () {return Edep;};
  void SetParticle (G4int i, G4double w, G4ThreeVector aparticle)
    {Slayer = i; Weight = w; Particle = aparticle;};

  G4ThreeVector GetParticle () {return Particle;};
  G4double GetWeight() {return Weight;};
  G4int GetLayer() {return Slayer;};

  void Draw () {};
  void Print () {};
  void EndOfEvent (G4HCofThisEvent* HCE) {};
  void clear () {};
  void DrawAll () {};
  void PrintAll () {};

private:
  G4double Edep;
  G4double Weight;
  G4ThreeVector Particle;
  G4int Slayer;

};

typedef G4THitsCollection<MLHit> MLHitsCollection;

extern G4Allocator<MLHit> MLHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* MLHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) MLHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void MLHit::operator delete(void *aHit)
{
  MLHitAllocator.FreeSingle((MLHit*) aHit);
}


////////////////////////////////////////////////////////////////////////////////
#endif

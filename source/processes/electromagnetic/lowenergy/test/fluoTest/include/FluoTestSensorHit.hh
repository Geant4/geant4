//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FluoTestSensorHit_h
#define FluoTestSensorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include <CLHEP/Random/Randomize.h>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FluoTestSensorHit : public G4VHit
{
public:

  FluoTestSensorHit();
  ~FluoTestSensorHit();
  FluoTestSensorHit(const FluoTestSensorHit&);
      const FluoTestSensorHit& operator=(const FluoTestSensorHit&);
      int operator==(const FluoTestSensorHit&) const;
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  void AddEnergy(G4double de)    {EdepTot += de;};

  void Draw();
  void Print();
  G4double GetEdepTot()      { return EdepTot;};
  G4double GetEdepDetect()   { return EdepDetect;};
  // G4double RandomCut (G4double energy);
  G4double RandomCut ();
  
private:

   G4double EdepTot;
  G4double Efficiency; 
  G4double EdepDetect; 
  G4double F;
  G4double deltaE;
  G4double epsilon;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4THitsCollection<FluoTestSensorHit> FluoTestSensorHitsCollection;

extern G4Allocator<FluoTestSensorHit> FluoTestSensorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* FluoTestSensorHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) FluoTestSensorHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void FluoTestSensorHit::operator delete(void* aHit)
{
  FluoTestSensorHitAllocator.FreeSingle((FluoTestSensorHit*) aHit);
}

#endif




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef fluoTestSensorHit_h
#define fluoTestSensorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class fluoTestSensorHit : public G4VHit
{
public:

  fluoTestSensorHit();
  ~fluoTestSensorHit();
  fluoTestSensorHit(const fluoTestSensorHit&);
  const fluoTestSensorHit& operator=(const fluoTestSensorHit&);
  int operator==(const fluoTestSensorHit&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);

  void Draw();
  void Print();
      
private:
  
  

public:
  
  void AddSi(G4double de, G4double dl) {EdepSi += de; TrackLengthSi += dl;};
  void AddSam(G4double de, G4double dl) {EdepSam += de; TrackLengthSam += dl;}; 
  void AddHPGe(G4double de, G4double dl) {EdepHPGe += de; TrackLengthHPGe += dl;};
 

  G4double GetEdepSi()     { return EdepSi; };
  G4double GetTrakSi()     { return TrackLengthSi; };
  G4double GetEdepHPGe()     { return EdepHPGe; };
  G4double GetTrakHPGe()     { return TrackLengthHPGe; };
  G4double GetEdepSam()     { return EdepSam; };
  G4double GetTrakSam()     { return TrackLengthSam; };

 
private:
  
  G4double EdepSi, TrackLengthSi;
  G4double EdepSam, TrackLengthSam;
  G4double EdepHPGe, TrackLengthHPGe;     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4THitsCollection<fluoTestSensorHit> fluoTestSensorHitsCollection;

extern G4Allocator<fluoTestSensorHit> fluoTestSensorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* fluoTestSensorHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) fluoTestSensorHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void fluoTestSensorHit::operator delete(void* aHit)
{
  fluoTestSensorHitAllocator.FreeSingle((fluoTestSensorHit*) aHit);
}

#endif




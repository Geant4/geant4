//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef XrayFluoSensorHit_h
#define XrayFluoSensorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoSensorHit : public G4VHit
{
public:
  
  XrayFluoSensorHit();
  ~XrayFluoSensorHit();
  XrayFluoSensorHit(const XrayFluoSensorHit&);
  const XrayFluoSensorHit& operator=(const XrayFluoSensorHit&);
  int operator==(const XrayFluoSensorHit&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  
  void Draw();
  void Print();
  
private:
  
  
  
public:
  
  void AddSen(G4double de, G4double dl) {EdepDet += de; TrackLengthDet += dl;};
  void AddSam(G4double de, G4double dl) {EdepSam += de; TrackLengthSam += dl;}; 
  
  
  G4double GetEdepDet()     { return EdepDet; };
  G4double GetTrakDet()     { return TrackLengthDet; };
  G4double GetEdepSam()     { return EdepSam; };
  G4double GetTrakSam()     { return TrackLengthSam; };
  
private:
  
  G4double EdepDet, TrackLengthDet;
  G4double EdepSam, TrackLengthSam;
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4THitsCollection<XrayFluoSensorHit> XrayFluoSensorHitsCollection;

extern G4Allocator<XrayFluoSensorHit> XrayFluoSensorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* XrayFluoSensorHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) XrayFluoSensorHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void XrayFluoSensorHit::operator delete(void* aHit)
{
  XrayFluoSensorHitAllocator.FreeSingle((XrayFluoSensorHit*) aHit);
}

#endif







// This Class describe the hits on the Payload

#ifndef GammaRayTelPayloadHit_h
#define GammaRayTelPayloadHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelPayloadHit : public G4VHit
{
public:
  
  GammaRayTelPayloadHit();
  ~GammaRayTelPayloadHit();
  GammaRayTelPayloadHit(const GammaRayTelPayloadHit&);
  const GammaRayTelPayloadHit& operator=(const GammaRayTelPayloadHit&);
  int operator==(const GammaRayTelPayloadHit&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  
  void Draw();
  void Print();

private:
  
  G4double EdepSil;  // Energy deposited on the silicon strip
  G4ThreeVector pos; // Position of the hit 
  G4int NStrip;      // Number of the strip
  G4int NSilPlane;   // Number of the plane
  G4int IsXPlane;    // Type of the plane (1 X, 0 Y)

public:
  
  inline void AddSil(G4double de) {EdepSil += de;};
  inline void SetNStrip(G4int i) {NStrip = i;};
  inline void SetNSilPlane(G4int i) {NSilPlane = i;};
  inline void SetPlaneType(G4int i) {IsXPlane = i;};
  inline void SetPos(G4ThreeVector xyz){ pos = xyz; }
  
  inline G4double GetEdepSil()     { return EdepSil; };
  inline G4int    GetNStrip()      { return NStrip; };
  inline G4int    GetNSilPlane()   { return NSilPlane; };
  inline G4int    GetPlaneType() {return IsXPlane;};      // Return 1 if the hit is on a X sil plane
  inline G4ThreeVector GetPos() { return pos; };
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4THitsCollection<GammaRayTelPayloadHit> GammaRayTelPayloadHitsCollection;

extern G4Allocator<GammaRayTelPayloadHit> GammaRayTelPayloadHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* GammaRayTelPayloadHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) GammaRayTelPayloadHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void GammaRayTelPayloadHit::operator delete(void* aHit)
{
  GammaRayTelPayloadHitAllocator.FreeSingle((GammaRayTelPayloadHit*) aHit);
}

#endif



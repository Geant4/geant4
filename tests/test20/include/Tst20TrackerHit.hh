// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20TrackerHit.hh,v 1.1 2001-05-24 19:49:21 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ Tst20TrackerHit  ------
// ************************************************************
// This Class describe the hits on the Tracker

#ifndef Tst20TrackerHit_h
#define Tst20TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst20TrackerHit : public G4VHit
{
public:
  
  Tst20TrackerHit();
  ~Tst20TrackerHit();
  Tst20TrackerHit(const Tst20TrackerHit&);
  const Tst20TrackerHit& operator=(const Tst20TrackerHit&);
  int operator==(const Tst20TrackerHit&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  
  void Draw();
  void Print();

private:
  
  G4double EdepSil;  // Energy deposited on the silicon strip
  G4ThreeVector pos; // Position of the hit 
  G4int NPixel;      // Number of the pixel
  G4int NSilPlane;   // Number of the plane

public:
  
  inline void AddSil(G4double de) {EdepSil += de;};
  inline void SetNPixel(G4int i) {NPixel = i;};
  inline void SetNSilPlane(G4int i) {NSilPlane = i;};
  inline void SetPos(G4ThreeVector xyz){ pos = xyz; }
  
  inline G4double GetEdepSil()     { return EdepSil; };
  inline G4int    GetNPixel()      { return NPixel; };
  inline G4int    GetNSilPlane()   { return NSilPlane; };
  inline G4ThreeVector GetPos() { return pos; };
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4THitsCollection<Tst20TrackerHit> Tst20TrackerHitsCollection;

extern G4Allocator<Tst20TrackerHit> Tst20TrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* Tst20TrackerHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) Tst20TrackerHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void Tst20TrackerHit::operator delete(void* aHit)
{
  Tst20TrackerHitAllocator.FreeSingle((Tst20TrackerHit*) aHit);
}

#endif










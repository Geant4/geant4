// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelTrackerHit.hh,v 1.1 2001-03-05 13:58:20 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ GammaRayTelTrackerHit  ------
//           by R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************
// This Class describe the hits on the Tracker

#ifndef GammaRayTelTrackerHit_h
#define GammaRayTelTrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelTrackerHit : public G4VHit
{
public:
  
  GammaRayTelTrackerHit();
  ~GammaRayTelTrackerHit();
  GammaRayTelTrackerHit(const GammaRayTelTrackerHit&);
  const GammaRayTelTrackerHit& operator=(const GammaRayTelTrackerHit&);
  int operator==(const GammaRayTelTrackerHit&) const;
  
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
  inline G4int    GetPlaneType()   {return IsXPlane;};      
  inline G4ThreeVector GetPos() { return pos; };
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4THitsCollection<GammaRayTelTrackerHit> GammaRayTelTrackerHitsCollection;

extern G4Allocator<GammaRayTelTrackerHit> GammaRayTelTrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* GammaRayTelTrackerHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) GammaRayTelTrackerHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void GammaRayTelTrackerHit::operator delete(void* aHit)
{
  GammaRayTelTrackerHitAllocator.FreeSingle((GammaRayTelTrackerHit*) aHit);
}

#endif










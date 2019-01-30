//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
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
  G4bool operator==(const GammaRayTelTrackerHit&) const;
  
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

extern G4ThreadLocal G4Allocator<GammaRayTelTrackerHit> *GammaRayTelTrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* GammaRayTelTrackerHit::operator new(size_t)
{
  if (!GammaRayTelTrackerHitAllocator)
    GammaRayTelTrackerHitAllocator = new G4Allocator<GammaRayTelTrackerHit>;
  return (void*) GammaRayTelTrackerHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void GammaRayTelTrackerHit::operator delete(void* aHit)
{
  GammaRayTelTrackerHitAllocator->FreeSingle((GammaRayTelTrackerHit*) aHit);
}

#endif










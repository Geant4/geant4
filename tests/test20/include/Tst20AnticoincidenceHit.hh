// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20AnticoincidenceHit.hh,v 1.1 2001-05-24 19:49:21 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ Tst20AnticoincidenceHit  ------
//
// ************************************************************
// This Class describe the hits on the Anticoincidence

#ifndef Tst20AnticoincidenceHit_h
#define Tst20AnticoincidenceHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst20AnticoincidenceHit : public G4VHit
{
public:
  
  Tst20AnticoincidenceHit();
  ~Tst20AnticoincidenceHit();
  Tst20AnticoincidenceHit(const Tst20AnticoincidenceHit&);
  const Tst20AnticoincidenceHit& operator=(const Tst20AnticoincidenceHit&);
  int operator==(const Tst20AnticoincidenceHit&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  
  void Draw();
  void Print();

private:
  
  G4double EdepACD;  // Energy deposited on the ACD tile
  G4ThreeVector pos; // Position of the hit
  G4int ACDTileNumber; // Number of the ACD tile
  G4int IsACDPlane;    // Type of the plane (0 top, 1 L-R, 2 F-R)

public:
  
  inline void AddEnergy(G4double de) {EdepACD += de;};
  inline void SetACDTileNumber(G4int i) {ACDTileNumber = i;};
  inline void SetACDType(G4int i) {IsACDPlane = i;};
  inline void SetPos(G4ThreeVector xyz){ pos = xyz; }
  
  inline G4double GetEdepACD()     { return EdepACD; };
  inline G4int    GetACDTileNumber()   { return ACDTileNumber; };
  inline G4int    GetACDType()   {return IsACDPlane;};      
  inline G4ThreeVector GetPos() { return pos; };
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4THitsCollection<Tst20AnticoincidenceHit> Tst20AnticoincidenceHitsCollection;

extern G4Allocator<Tst20AnticoincidenceHit> Tst20AnticoincidenceHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* Tst20AnticoincidenceHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) Tst20AnticoincidenceHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void Tst20AnticoincidenceHit::operator delete(void* aHit)
{
  Tst20AnticoincidenceHitAllocator.FreeSingle((Tst20AnticoincidenceHit*) aHit);
}

#endif











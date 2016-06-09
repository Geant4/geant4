//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// --------------------------------------------------------------
//                 GEANT 4 - ULTRA experiment example
// --------------------------------------------------------------
//
// Code developed by:
// B. Tome, M.C. Espirito-Santo, A. Trindade, P. Rodrigues 
//
//   **********************************************
//   *        UltraOpticalHit.hh
//   **********************************************
//    
//    Classe defining an optical photon hit in the photomultiplier
//    Data members are the photon energy and incidence position 
//
#ifndef UltraOpticalHit_h
#define UltraOpticalHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class UltraOpticalHit : public G4VHit
{
public:
  
  UltraOpticalHit();
  ~UltraOpticalHit();
  UltraOpticalHit(const UltraOpticalHit&);
  const UltraOpticalHit& operator=(const UltraOpticalHit&);
  int operator==(const UltraOpticalHit&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  
  void Draw();
  void Print();

private:
  
  G4double      fPhotEne;            // Photon energy energy 
  G4ThreeVector fPhotPos;            // Position of the hit

public:

  // Set functions to store information on hits  
  inline void SetEnergy(G4double fEn)   {fPhotEne= fEn;}
  inline void SetPosition(G4ThreeVector xyz) {fPhotPos = xyz;}

  // Get functions to acess information on hits
  inline G4double      GetEnergy()    {return fPhotEne; }
  inline G4ThreeVector GetPosition()  {return fPhotPos; } 

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4THitsCollection<UltraOpticalHit> UltraOpticalHitsCollection;

extern G4Allocator<UltraOpticalHit> UltraOpticalHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* UltraOpticalHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) UltraOpticalHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void UltraOpticalHit::operator delete(void* aHit)
{
  UltraOpticalHitAllocator.FreeSingle((UltraOpticalHit*) aHit);
}

#endif

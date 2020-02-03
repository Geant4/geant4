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
  G4bool operator==(const UltraOpticalHit&) const;
  
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

extern G4ThreadLocal G4Allocator<UltraOpticalHit> *UltraOpticalHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* UltraOpticalHit::operator new(size_t)
{
  if (!UltraOpticalHitAllocator)
    UltraOpticalHitAllocator = new G4Allocator<UltraOpticalHit>;
  return (void*) UltraOpticalHitAllocator->MallocSingle();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void UltraOpticalHit::operator delete(void* aHit)
{
  UltraOpticalHitAllocator->FreeSingle((UltraOpticalHit*) aHit);
}

#endif

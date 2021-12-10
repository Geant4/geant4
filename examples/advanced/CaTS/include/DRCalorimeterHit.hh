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

// ********************************************************************
//
//  CaTS (Calorimetry and Tracking Simulation)
//
//  Authors : Hans Wenzel
//            Soon Yung Jun
//            (Fermi National Accelerator Laboratory)
//
// History
//   October 18th, 2021 : first implementation
//
// ********************************************************************
//
/// \file DRCalorimeterHit.hh
/// \brief Definition of the CaTS::DRCalorimeterHit class

#pragma once

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class DRCalorimeterHit : public G4VHit
{
 public:
  DRCalorimeterHit();
  ~DRCalorimeterHit() = default;
  DRCalorimeterHit(const DRCalorimeterHit&);
  const DRCalorimeterHit& operator=(const DRCalorimeterHit&);
  G4bool operator==(const DRCalorimeterHit&) const;
  inline void* operator new(size_t);
  inline void operator delete(void*);
  void Draw() final;
  inline void Print() final
  {
    G4cout << "DRCalorimeterHit  id:  " << fid << " Edep: " << fEdep
           << " em_Edep: " << fem_Edep << " NCeren: " << fNceren
           << " X: " << fposition.getX() << " Y: " << fposition.getY()
           << " Z: " << fposition.getZ() << G4endl;
  }

  DRCalorimeterHit(unsigned int i, G4double e, G4double em, unsigned int nc,
                   G4double t, G4ThreeVector p);
  inline void SetPosition(G4ThreeVector position) { fposition = position; };
  inline G4ThreeVector GetPosition() const { return fposition; };
  inline void SetTime(G4double time) { ftime = time; };
  inline G4double GetTime() const { return ftime; };
  inline void SetId(unsigned int id) { fid = id; };
  inline unsigned int GetId() const { return fid; };
  inline void SetEdep(G4double Edep) { fEdep = Edep; };
  inline G4double GetEdep() const { return fEdep; };
  inline void SetNceren(unsigned int Nceren) { fNceren = Nceren; };
  inline unsigned int GetNceren() const { return fNceren; };
  inline void SetEm_Edep(G4double em_Edep) { fem_Edep = em_Edep; };
  inline G4double GetEm_Edep() const { return fem_Edep; };

 private:
  unsigned int fid{ 0 };
  G4double fEdep{ 0 };
  G4double fem_Edep{ 0 };
  unsigned int fNceren{ 0 };
  G4double ftime{ 0 };
  G4ThreeVector fposition{ 0, 0, 0 };
};

using DRCalorimeterHitsCollection = G4THitsCollection<DRCalorimeterHit>;
extern G4ThreadLocal G4Allocator<DRCalorimeterHit>* DRCalorimeterHitAllocator;

inline void* DRCalorimeterHit::operator new(size_t)
{
  if(!DRCalorimeterHitAllocator)
  {
    DRCalorimeterHitAllocator = new G4Allocator<DRCalorimeterHit>;
  }
  return (void*) DRCalorimeterHitAllocator->MallocSingle();
}

inline void DRCalorimeterHit::operator delete(void* aHit)
{
  DRCalorimeterHitAllocator->FreeSingle((DRCalorimeterHit*) aHit);
}

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
/// \file CalorimeterHit.hh
/// \brief Definition of the CaTS::CalorimeterHit class

#pragma once

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class CalorimeterHit : public G4VHit
{
 public:
  CalorimeterHit();
  ~CalorimeterHit() = default;
  CalorimeterHit(const CalorimeterHit&);
  const CalorimeterHit& operator=(const CalorimeterHit&);
  G4bool operator==(const CalorimeterHit&) const;
  inline void* operator new(size_t);
  inline void operator delete(void*);
  void Draw() final;
  inline void Print() final
  {
    G4cout << "CalorimeterHit  id:  " << fid << " Edep: " << fEdep
           << " em_Edep: " << fem_Edep << " time: " << ftime
           << " X: " << fposition.getX() << " Y: " << fposition.getY()
           << " Z: " << fposition.getZ() << G4endl;
  }
  CalorimeterHit(unsigned int i, G4double e, G4double em, G4double t,
                 G4ThreeVector p);
  inline void SetPosition(G4ThreeVector position) { fposition = position; };
  inline G4ThreeVector GetPosition() const { return fposition; };
  inline void SetTime(G4double time) { ftime = time; };
  inline G4double GetTime() const { return ftime; };
  inline void SetId(unsigned int id) { fid = id; };
  inline unsigned int GetId() const { return fid; };
  inline void SetEdep(G4double Edep) { fEdep = Edep; };
  inline G4double GetEdep() const { return fEdep; };
  inline void Setem_Edep(G4double em_Edep) { fem_Edep = em_Edep; };
  inline G4double Getem_Edep() const { return fem_Edep; };

 private:
  unsigned int fid{ 0 };
  G4double fEdep{ 0 };
  G4double fem_Edep{ 0 };
  G4double ftime{ 0 };
  G4ThreeVector fposition{ 0, 0, 0 };
};

using CalorimeterHitsCollection = G4THitsCollection<CalorimeterHit>;
extern G4ThreadLocal G4Allocator<CalorimeterHit>* CalorimeterHitAllocator;

inline void* CalorimeterHit::operator new(size_t)
{
  if(!CalorimeterHitAllocator)
  {
    CalorimeterHitAllocator = new G4Allocator<CalorimeterHit>;
  }
  return (void*) CalorimeterHitAllocator->MallocSingle();
}

inline void CalorimeterHit::operator delete(void* aHit)
{
  CalorimeterHitAllocator->FreeSingle((CalorimeterHit*) aHit);
}

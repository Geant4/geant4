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
/// \file MscHit.hh
/// \brief Definition of the CaTS::MscHit class

#pragma once

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include <G4Types.hh>
#include <G4ios.hh>
#include <ostream>
#include <tls.hh>

class MscHit : public G4VHit
{
 public:
  MscHit();
  ~MscHit() = default;
  MscHit(const MscHit&);
  const MscHit& operator=(const MscHit&);
  G4bool operator==(const MscHit&) const;
  inline void* operator new(size_t);
  inline void operator delete(void*);
  void Draw() final;
  inline void Print() final
  {
    G4cout << "MscHit KinE: " << fkinE << " Px: " << fmomentum.getX()
           << " Py: " << fmomentum.getY() << " Pz: " << fmomentum.getZ()
           << G4endl;
  }
  MscHit(G4double kinE, G4ThreeVector momentum);
  inline void SetKinE(G4double kinE) { fkinE = kinE; };
  inline G4double GetKinE() const { return fkinE; };
  inline void SetMomentum(G4ThreeVector momentum) { fmomentum = momentum; };
  inline G4ThreeVector GetMomentum() const { return fmomentum; };

 private:
  G4double fkinE{ 0 };
  G4ThreeVector fmomentum{ 0, 0, 0 };
};

using MscHitsCollection = G4THitsCollection<MscHit>;
extern G4ThreadLocal G4Allocator<MscHit>* MscHitAllocator;

inline void* MscHit::operator new(size_t)
{
  if(!MscHitAllocator)
  {
    MscHitAllocator = new G4Allocator<MscHit>;
  }
  return (void*) MscHitAllocator->MallocSingle();
}

inline void MscHit::operator delete(void* aHit)
{
  MscHitAllocator->FreeSingle((MscHit*) aHit);
}

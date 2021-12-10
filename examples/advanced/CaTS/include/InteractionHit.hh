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
/// \file InteractionHit.hh
/// \brief Definition of the CaTS::InteractionHit class

#pragma once

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

class InteractionHit : public G4VHit
{
 private:
  G4String fpname{ "" };  // name of secondary particle
  G4double fpmom{ 0 };    // momentum of secondary particle
  G4double fEkin{ 0 };    // kinetic energy of secondary particle
  G4double ftheta{ 0 };   // theta of secondary particle
 public:
  InteractionHit();
  InteractionHit(G4String n, G4double m, G4double e, G4double t);
  ~InteractionHit() = default;
  InteractionHit(const InteractionHit&);
  const InteractionHit& operator=(const InteractionHit&);
  G4int operator==(const InteractionHit&) const;
  inline void* operator new(size_t);
  inline void operator delete(void*);
  inline void Print() final
  {
    G4cout << "InteractionHit pname : " << fpname
           << "  momentum [GeV]: " << fpmom << "  kinetic Energy [GeV]" << fEkin
           << "  theta: " << ftheta << G4endl;
  }
  inline void SetPname(G4String de) { fpname = de; };
  inline void SetPmom(G4double de) { fpmom = de; };
  inline void SetEkin(G4double de) { fEkin = de; };
  inline void SetTheta(G4double de) { ftheta = de; };
  inline G4String GetPname() { return fpname; };
  inline G4double GetPmom() { return fpmom; };
  inline G4double GetEkin() { return fEkin; };
  inline G4double GetTheta() { return ftheta; };
};

using InteractionHitsCollection = G4THitsCollection<InteractionHit>;
extern G4ThreadLocal G4Allocator<InteractionHit>* InteractionHitAllocator;

inline void* InteractionHit::operator new(size_t)
{
  if(!InteractionHitAllocator)
  {
    InteractionHitAllocator = new G4Allocator<InteractionHit>;
  }
  return (void*) InteractionHitAllocator->MallocSingle();
}

inline void InteractionHit::operator delete(void* aHit)
{
  InteractionHitAllocator->FreeSingle((InteractionHit*) aHit);
}

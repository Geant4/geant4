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
/// \file ExGflashHit.hh
/// \brief Definition of the ExGflashHit class
//
#ifndef ExGflashHit_h
#define ExGflashHit_h 1

#include "G4Allocator.hh"
#include "G4RotationMatrix.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"

class G4LogicalVolume;

class ExGflashHit : public G4VHit
{
 public:
  ExGflashHit();
  ~ExGflashHit() override;
  ExGflashHit(const ExGflashHit& right);
  ExGflashHit& operator=(const ExGflashHit& right);
  G4bool operator==(const ExGflashHit& right) const;

  inline void* operator new(size_t);
  inline void operator delete(void* aHit);
  void* operator new(size_t, void* p) { return p; }
#ifndef G4NOT_ISO_DELETES
  void operator delete(void*, void*) {}
#endif

  void Draw() override;
  void Print() override;

 private:
  G4double fEdep;
  G4ThreeVector fPos;

 public:
  void SetEdep(G4double de)
  {
    fEdep = de;
  };
  void AddEdep(G4double de)
  {
    fEdep += de;
  };
  G4double GetEdep()
  {
    return fEdep;
  };
  void SetPos(G4ThreeVector xyz)
  {
    fPos = xyz;
  };
  G4ThreeVector GetPos()
  {
    return fPos;
  };
};

using ExGflashHitsCollection = G4THitsCollection<ExGflashHit>;

extern G4ThreadLocal G4Allocator<ExGflashHit>* ExGflashHitAllocator;

inline void* ExGflashHit::operator new(size_t)
{
  if ( ExGflashHitAllocator == nullptr )
    ExGflashHitAllocator = new G4Allocator<ExGflashHit>;
  return (void*)ExGflashHitAllocator->MallocSingle();
}

inline void ExGflashHit::operator delete(void* aHit)
{
  ExGflashHitAllocator->FreeSingle((ExGflashHit*)aHit);
}

#endif

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
/// \file RE05/include/RE05MuonHit.hh
/// \brief Definition of the RE05MuonHit class
//

#ifndef RE05MuonHit_h
#define RE05MuonHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"

class G4AttDef;

class RE05MuonHit : public G4VHit
{
  public:

      RE05MuonHit();
      virtual ~RE05MuonHit();
      RE05MuonHit(const RE05MuonHit &right);
      const RE05MuonHit& operator=(const RE05MuonHit &right);
      G4bool operator==(const RE05MuonHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      virtual void Draw();
      virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
      virtual std::vector<G4AttValue>* CreateAttValues() const;
      virtual void Print();

  private:
      G4double fEdep;
      G4ThreeVector fPos;
      static std::map<G4String,G4AttDef> fAttDefs;

  public:
      inline void SetEdep(G4double de)
      { fEdep = de; }
      inline void AddEdep(G4double de)
      { fEdep += de; }
      inline G4double GetEdep()
      { return fEdep; }
      inline void SetPos(G4ThreeVector xyz)
      { fPos = xyz; }
      inline G4ThreeVector GetPos()
      { return fPos; }
};

typedef G4THitsCollection<RE05MuonHit> RE05MuonHitsCollection;

extern G4ThreadLocal G4Allocator<RE05MuonHit>* RE05MuonHitAllocator;

inline void* RE05MuonHit::operator new(size_t)
{
  if(!RE05MuonHitAllocator) RE05MuonHitAllocator = new G4Allocator<RE05MuonHit>;
  return (void *) RE05MuonHitAllocator->MallocSingle();
}

inline void RE05MuonHit::operator delete(void *aHit)
{
  RE05MuonHitAllocator->FreeSingle((RE05MuonHit*) aHit);
}

#endif

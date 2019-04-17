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
/// \file RE05/include/RE05TrackerHit.hh
/// \brief Definition of the RE05TrackerHit class
//

#ifndef RE05TrackerHit_h
#define RE05TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"

class G4AttDef;

class RE05TrackerHit : public G4VHit
{
  public:

      RE05TrackerHit();
      virtual ~RE05TrackerHit();
      RE05TrackerHit(const RE05TrackerHit &right);
      const RE05TrackerHit& operator=(const RE05TrackerHit &right);
      G4bool operator==(const RE05TrackerHit &right) const;

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
      inline G4double GetEdep()
      { return fEdep; }
      inline void SetPos(G4ThreeVector xyz)
      { fPos = xyz; }
      inline G4ThreeVector GetPos()
      { return fPos; }
};

typedef G4THitsCollection<RE05TrackerHit> RE05TrackerHitsCollection;

extern G4ThreadLocal G4Allocator<RE05TrackerHit>* RE05TrackerHitAllocator;

inline void* RE05TrackerHit::operator new(size_t)
{
  if(!RE05TrackerHitAllocator) RE05TrackerHitAllocator = new G4Allocator<RE05TrackerHit>;
  return (void *) RE05TrackerHitAllocator->MallocSingle();
}

inline void RE05TrackerHit::operator delete(void *aHit)
{
  RE05TrackerHitAllocator->FreeSingle((RE05TrackerHit*) aHit);
}

#endif

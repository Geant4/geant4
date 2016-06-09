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
// $Id: RE01TrackerHit.hh,v 1.2 2005/06/02 21:30:50 perl Exp $
// GEANT4 tag $Name: geant4-08-00 $
//


#ifndef RE01TrackerHit_h
#define RE01TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class G4AttDef;
class G4AttValue;

class RE01TrackerHit : public G4VHit
{
  public:

      RE01TrackerHit();
      ~RE01TrackerHit();
      RE01TrackerHit(const RE01TrackerHit &right);
      const RE01TrackerHit& operator=(const RE01TrackerHit &right);
      G4int operator==(const RE01TrackerHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      void Draw();
      virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
      virtual std::vector<G4AttValue>* CreateAttValues() const;
      void Print();

  private:
      G4double edep;
      G4ThreeVector pos;
      G4int trackID;

  public:
      inline void SetEdep(G4double de)
      { edep = de; }
      inline G4double GetEdep() const
      { return edep; }
      inline void SetPos(G4ThreeVector xyz)
      { pos = xyz; }
      inline G4ThreeVector GetPos() const
      { return pos; }
      inline void SetTrackID(G4int i)
      { trackID = i; }
      inline G4int GetTrackID() const
      { return trackID; }
};

typedef G4THitsCollection<RE01TrackerHit> RE01TrackerHitsCollection;

extern G4Allocator<RE01TrackerHit> RE01TrackerHitAllocator;

inline void* RE01TrackerHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) RE01TrackerHitAllocator.MallocSingle();
  return aHit;
}

inline void RE01TrackerHit::operator delete(void *aHit)
{
  RE01TrackerHitAllocator.FreeSingle((RE01TrackerHit*) aHit);
}

#endif

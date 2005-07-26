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

#ifndef G4HumanPhantomHit_h
#define G4HumanPhantomHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"

class G4HumanPhantomHit : public G4VHit
{
  public:

      G4HumanPhantomHit();
     ~G4HumanPhantomHit();
      G4HumanPhantomHit(const G4HumanPhantomHit&);
      const G4HumanPhantomHit& operator=(const G4HumanPhantomHit&);
      G4int operator==(const G4HumanPhantomHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Draw();
      void Print();

  public:
  
      void SetTrackID  (G4int track)      { trackID = track; };
      void SetBodyPartID(G4int bpID)      { bodypartID = bpID; };

      void SetBodyPartName(G4String bpName)      { bodypartName = bpName; };

      void SetEdep     (G4double de)      { edep = de; };
      void SetPos      (G4ThreeVector xyz){ pos = xyz; };
      
      G4int GetTrackID()    { return trackID; };
      G4int GetBodyPartID() { return bodypartID; };

      G4String GetBodyPartName() { return bodypartName; };

      G4double GetEdep()    { return edep; };      
      G4ThreeVector GetPos(){ return pos; };
      
  private:
  
      G4int         trackID;
      G4int         bodypartID;

      G4String      bodypartName;

      G4double      edep;
      G4ThreeVector pos;
};


typedef G4THitsCollection<G4HumanPhantomHit> G4HumanPhantomHitsCollection;

extern G4Allocator<G4HumanPhantomHit> G4HumanPhantomHitAllocator;


inline void* G4HumanPhantomHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) G4HumanPhantomHitAllocator.MallocSingle();
  return aHit;
}


inline void G4HumanPhantomHit::operator delete(void *aHit)
{
  G4HumanPhantomHitAllocator.FreeSingle((G4HumanPhantomHit*) aHit);
}


#endif

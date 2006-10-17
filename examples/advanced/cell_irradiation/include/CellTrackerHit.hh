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
//    **************************************
//    *                                    *
//    *        CellTRackerHit.hh           *
//    *                                    *
//    **************************************
//
//

#ifndef CellTrackerHit_h
#define CellTrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class CellTrackerHit : public G4VHit
{
  public:

      CellTrackerHit();
     ~CellTrackerHit();
      CellTrackerHit(const CellTrackerHit&);
      const CellTrackerHit& operator=(const CellTrackerHit&);
      G4int operator==(const CellTrackerHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Draw();
      void Print();

  public:
  // Store the energy deposit and position in the hit
      void SetEdep    (G4double de)      { edep = de; };
      void SetPos     (G4int voxel){ voxel_hit = voxel; };
     
      G4double GetEdep()    { return edep; };      
      G4int GetPos(){ return voxel_hit; };
      
  private:
      G4double      edep;
      G4int voxel_hit;
};

typedef G4THitsCollection<CellTrackerHit> CellTrackerHitsCollection;

extern G4Allocator<CellTrackerHit> CellTrackerHitAllocator;

inline void* CellTrackerHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) CellTrackerHitAllocator.MallocSingle();
  return aHit;
}

inline void CellTrackerHit::operator delete(void *aHit)
{
  CellTrackerHitAllocator.FreeSingle((CellTrackerHit*) aHit);
}

#endif

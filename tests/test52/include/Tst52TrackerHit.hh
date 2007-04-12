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
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// $Id: Tst52TrackerHit.hh,v 1.1 2007-04-12 12:00:17 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst52TrackerHit_h
#define Tst52TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class Tst52TrackerHit : public G4VHit
{
  public:

      Tst52TrackerHit();
     ~Tst52TrackerHit();
      Tst52TrackerHit(const Tst52TrackerHit&);
      const Tst52TrackerHit& operator=(const Tst52TrackerHit&);
      G4int operator==(const Tst52TrackerHit&) const;

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

typedef G4THitsCollection<Tst52TrackerHit> Tst52TrackerHitsCollection;

extern G4Allocator<Tst52TrackerHit> Tst52TrackerHitAllocator;

inline void* Tst52TrackerHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) Tst52TrackerHitAllocator.MallocSingle();
  return aHit;
}

inline void Tst52TrackerHit::operator delete(void *aHit)
{
  Tst52TrackerHitAllocator.FreeSingle((Tst52TrackerHit*) aHit);
}

#endif

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

#ifndef ExN04MuonHit_h
#define ExN04MuonHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class ExN04MuonHit : public G4VHit
{
  public:

      ExN04MuonHit();
      ~ExN04MuonHit();
      ExN04MuonHit(const ExN04MuonHit &right);
      const ExN04MuonHit& operator=(const ExN04MuonHit &right);
      int operator==(const ExN04MuonHit &right) const;


      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      void Draw();
      void Print();

  private:
      G4double edep;
      G4ThreeVector pos;

  public:
      inline void SetEdep(G4double de)
      { edep = de; }
      inline void AddEdep(G4double de)
      { edep += de; }
      inline G4double GetEdep()
      { return edep; }
      inline void SetPos(G4ThreeVector xyz)
      { pos = xyz; }
      inline G4ThreeVector GetPos()
      { return pos; }

};

typedef G4THitsCollection<ExN04MuonHit> ExN04MuonHitsCollection;

extern G4Allocator<ExN04MuonHit> ExN04MuonHitAllocator;

inline void* ExN04MuonHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) ExN04MuonHitAllocator.MallocSingle();
  return aHit;
}

inline void ExN04MuonHit::operator delete(void *aHit)
{
  ExN04MuonHitAllocator.FreeSingle((ExN04MuonHit*) aHit);
}

#endif



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

#ifndef RemSimHit_h
#define RemSimHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class RemSimHit : public G4VHit
{
public:

  RemSimHit();
  ~RemSimHit();
  RemSimHit(const RemSimHit &right);
  const RemSimHit& operator=(const RemSimHit &right);
  G4int operator==(const RemSimHit &right) const;

  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

  void Draw();
  void Print();

private:
  G4double edep;
  G4int xID;
  G4int yID;
  G4int zID;
  G4ThreeVector position;

public:
  inline void SetEdep(G4double de)
  { edep = de; }
  inline G4double GetEdep()
  { return edep; }
  inline  void SetIndexes(G4int i)
  {zID = i;}
 
  inline  G4int GetIndexZ(){return zID;} 
  
  inline void SetPosition(G4ThreeVector xyz)
      { position = xyz; }
      
  inline G4ThreeVector GetPosition(){ return position; }
};

typedef G4THitsCollection<RemSimHit> RemSimHitsCollection;

extern G4Allocator<RemSimHit> RemSimHitAllocator;

inline void* RemSimHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) RemSimHitAllocator.MallocSingle();
  return aHit;
}

inline void RemSimHit::operator delete(void *aHit)
{
  RemSimHitAllocator.FreeSingle((RemSimHit*) aHit);
}

#endif

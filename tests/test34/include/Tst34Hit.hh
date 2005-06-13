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
#ifndef Tst34Hit_h
#define Tst34Hit_h 1
 
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class Tst34Hit : public G4VHit
{
  public:

      Tst34Hit();
      Tst34Hit(G4LogicalVolume* logVol);
      ~Tst34Hit();
      Tst34Hit(const Tst34Hit &right);
      const Tst34Hit& operator=(const Tst34Hit &right);
      int operator==(const Tst34Hit &right) const;


      inline void *operator new(size_t);
      inline void operator delete(void *aHit);
      void *operator new(size_t,void*p){return p;}
#ifndef G4NOT_ISO_DELETES
      void operator delete(void*,void*){}
#endif

      void Draw();
      void Print();

  private:

      G4double edep;
      G4ThreeVector pos;
      G4int crystalnumber;
      G4ThreeVector start;
      G4RotationMatrix rot;
      const G4LogicalVolume* pLogV;

  public:

      inline void SetEdep(G4double de)
      { edep = de; }
      inline void AddEdep(G4double de)
      { edep += de; }
      inline G4double GetEdep()
      { return edep; }
      inline void SetPos(G4ThreeVector xyz)
      { pos = xyz; }
      inline G4int GetCrystalNum()
      { return crystalnumber; }
      inline void SetCrystalNum(G4int num)
      { crystalnumber=num; }
      inline G4ThreeVector GetPos()
      { return pos; }
      inline void SetStart(G4ThreeVector xyz)
      { start = xyz; }
      inline G4ThreeVector GetStart()
      { return start; }

      inline void SetRot(G4RotationMatrix rmat)
      { rot = rmat; }
      inline G4RotationMatrix GetRot()
      { return rot; }
      inline const G4LogicalVolume * GetLogV()
      { return pLogV; }
};

typedef G4THitsCollection<Tst34Hit> Tst34HitsCollection;

extern G4Allocator<Tst34Hit> Tst34HitAllocator;

inline void* Tst34Hit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) Tst34HitAllocator.MallocSingle();
  return aHit;
}

inline void Tst34Hit::operator delete(void *aHit)
{
  Tst34HitAllocator.FreeSingle((Tst34Hit*) aHit);
}

#endif

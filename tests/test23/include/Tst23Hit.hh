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
//
// $Id: Tst23Hit.hh,v 1.1 2001-12-14 14:53:23 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst23Hit_h
#define Tst23Hit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class Tst23Hit : public G4VHit
{
  public:

      Tst23Hit();
      Tst23Hit(G4int id);
      ~Tst23Hit();
      Tst23Hit(const Tst23Hit &right);
      const Tst23Hit& operator=(const Tst23Hit &right);
      int operator==(const Tst23Hit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);
      void *operator new(size_t,void*p){return p;}
#ifndef G4NOT_ISO_DELETES
      void operator delete(void *aHit,void*){}
#endif

      void Draw();
      void Print();

  private:
      G4double edep;
      G4int id;

  public:
      inline void SetEdep(G4double de)
      { edep = de; };
      inline void AddEdep(G4double de)
      { edep += de; };
      inline G4double GetEdep()
      { return edep; };
     inline G4double GetID()
      { return id; };

};

typedef G4THitsCollection<Tst23Hit> Tst23HitsCollection;

extern G4Allocator<Tst23Hit> Tst23HitAllocator;

inline void* Tst23Hit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) Tst23HitAllocator.MallocSingle();
  return aHit;
}

inline void Tst23Hit::operator delete(void *aHit)
{
  Tst23HitAllocator.FreeSingle((Tst23Hit*) aHit);
}

#endif




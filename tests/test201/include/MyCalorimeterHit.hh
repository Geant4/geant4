// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyCalorimeterHit.hh,v 1.3 2000-05-26 13:11:39 barrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef MyCalorimeterHit_h
#define MyCalorimeterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class MyCalorimeterHit : public G4VHit
{
  public:

      MyCalorimeterHit();
      ~MyCalorimeterHit();
      MyCalorimeterHit(const MyCalorimeterHit &right);
      MyCalorimeterHit(G4VPhysicalVolume* physVol);
      const MyCalorimeterHit& operator=(const MyCalorimeterHit &right);
      int operator==(const MyCalorimeterHit &right) const;


      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      void Draw();
      void Print();

  private:
      G4double edep;
      G4ThreeVector pos;
      G4RotationMatrix rot;
      G4VPhysicalVolume* pPhys;

  public:
      inline void SetEdep(G4double de)
      { edep = de; };
      inline void AddEdep(G4double de)
      { edep += de; };
      inline G4double GetEdep()
      { return edep; };
      inline void SetPos(G4ThreeVector xyz)
      { pos = xyz; };
      inline G4ThreeVector GetPos()
      { return pos; };
      inline void SetRot(G4RotationMatrix rmat)
      { rot = rmat; };
      inline G4RotationMatrix GetRot()
      { return rot; };
      inline G4VPhysicalVolume* GetPPhys()
      { return pPhys; };

};

typedef G4THitsCollection<MyCalorimeterHit> MyCalorimeterHitsCollection;

extern G4Allocator<MyCalorimeterHit> MyCalorimeterHitAllocator;

inline void* MyCalorimeterHit::operator new(size_t)
{
  return (void *) MyCalorimeterHitAllocator.MallocSingle();
}

inline void MyCalorimeterHit::operator delete(void *aHit)
{
  MyCalorimeterHitAllocator.FreeSingle((MyCalorimeterHit*) aHit);
}

#endif



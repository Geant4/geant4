// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LogicalVolumeStore.hh,v 1.2 1999-11-11 15:35:49 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4LogicalVolumeStore
//
// Container for all LogicalVolumes, with functionality derived from
// G4RWTPtrOrderedVector<T>. The class is `singleton', in that only
// one can exist, and access is facillitated via G4LogicalVolumeStore::GetInstance()
//
// All LogicalVolumes should be registered with G4LogicalVolumeStore, and removed on their
// destruction. Intended principally for UI browser. The underlying
// container initially has a capacity of 100.
//
// If much additional functionality is added, should consider containment
// instead of inheritance for G4RWTPtrOrderedVector<T>
//
// Class member functions:
//
// static void Register(G4LogicalVolume* pVolume)
//   Add the logical volume to the collection
// static void DeRegister(G4LogicalVolume* pVolume)
//   REmove the logical volume from the collection
// static G4LogicalVolumeStore* GetInstance()
//   Get a ptr to the unique G4LogicalVolumeStore, creting it if necessary
//
// Member functions:
//
// [as per RWTPtrOrderedvector]
//
// NOTE: Constructor is protected - creation and subsequent access is via
//       GetInstance
//
// Member data:
//
// static G4LogicalVolumeStore* fgInstance
//   Ptr to the single G4LogicalVolumeStore
//
// History:
// 10.07.95 P.Kent Initial version

#ifndef G4VLOGICALVOLUMESTORE_HH
#define G4VLOGICALVOLUMESTORE_HH

#include "g4rw/tpordvec.h"

#include "G4LogicalVolume.hh"

class G4LogicalVolumeStore : public G4RWTPtrOrderedVector<G4LogicalVolume>
{
  public:
    static void Register(G4LogicalVolume* pVolume);
    static void DeRegister(G4LogicalVolume* pVolume);
    static G4LogicalVolumeStore* GetInstance();
    virtual ~G4LogicalVolumeStore();
  protected:
    G4LogicalVolumeStore();
  private:
    static G4LogicalVolumeStore* fgInstance;
};

#endif

// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LogicalVolumeStore.hh,v 1.4 2000-04-20 16:49:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4LogicalVolumeStore
//
// Class description:
//
// Container for all LogicalVolumes, with functionality derived from
// G4RWTPtrOrderedVector<T>. The class is a `singleton', in that only
// one can exist, and access is provided via the static function
// G4LogicalVolumeStore::GetInstance()
//
// All LogicalVolumes should be registered with G4LogicalVolumeStore,
// and removed on their destruction. Intended principally for UI browser.
// The underlying container initially has a capacity of 100.
//
// If much additional functionality is added, should consider containment
// instead of inheritance for G4RWTPtrOrderedVector<T>.
//
// Member data:
//
// static G4LogicalVolumeStore* fgInstance
//   - Ptr to the single G4LogicalVolumeStore.

// History:
// 10.07.95 P.Kent Initial version

#ifndef G4VLOGICALVOLUMESTORE_HH
#define G4VLOGICALVOLUMESTORE_HH

#include "g4rw/tpordvec.h"

#include "G4LogicalVolume.hh"

class G4LogicalVolumeStore : public G4RWTPtrOrderedVector<G4LogicalVolume>
{
  public:  // with description

    static void Register(G4LogicalVolume* pVolume);
      // Add the logical volume to the collection.
    static void DeRegister(G4LogicalVolume* pVolume);
      // Remove the logical volume from the collection.
    static G4LogicalVolumeStore* GetInstance();
      // Get a ptr to the unique G4LogicalVolumeStore,
      // creating it if necessary.

    virtual ~G4LogicalVolumeStore();
      // Destructor.

  protected:

    G4LogicalVolumeStore();

  private:

    static G4LogicalVolumeStore* fgInstance;
};

#endif

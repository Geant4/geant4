// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicalVolumeStore.hh,v 1.6 2001-04-20 20:13:53 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4PhysicalVolume
//
// Class description:
//
// Container for all solids, with functionality derived from
// std::vector<T>. The class is a `singleton', in that only
// one can exist, and access is provided via the static method
// G4PhysicalVolumeStore::GetInstance()
//
// All solids should be registered with G4PhysicalVolumeStore, and removed on
// their destruction. Intended principally for UI browser. The underlying
// container initially has a capacity of 100.
//
// If much additional functionality is added, should consider containment
// instead of inheritance for std::vector<T>
//
// Member data:
//
// static G4PhysicalVolumeStore*
//   - Ptr to the single G4PhysicalVolumeStore.

// History:
// 18.04.01 G.Cosmo Migrated to STL vector
// 25.07.95 P.Kent  Initial version

#ifndef G4PHYSICALVOLUMESTORE_HH
#define G4PHYSICALVOLUMESTORE_HH

#include "g4std/vector"

#include "G4VPhysicalVolume.hh"

class G4PhysicalVolumeStore : public G4std::vector<G4VPhysicalVolume*>
{
  public:  // with description

    static void Register(G4VPhysicalVolume* pSolid);
      // Add the volume to the collection.
    static void DeRegister(G4VPhysicalVolume* pSolid);
      // Remove the volume from the collection.
    static G4PhysicalVolumeStore* GetInstance();
      // Get a ptr to the unique G4PhysicalVolumeStore,
      // creating it if necessary.

    virtual ~G4PhysicalVolumeStore();
      // Default destructor.

  protected:

    G4PhysicalVolumeStore();

  private:

    static G4PhysicalVolumeStore* fgInstance;
};

#endif

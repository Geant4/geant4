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
// $Id: G4LogicalVolumeStore.hh,v 1.6 2001-07-11 09:59:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4LogicalVolumeStore
//
// Class description:
//
// Container for all LogicalVolumes, with functionality derived from
// std::vector<T>. The class is a `singleton', in that only
// one can exist, and access is provided via the static function
// G4LogicalVolumeStore::GetInstance()
//
// All LogicalVolumes should be registered with G4LogicalVolumeStore,
// and removed on their destruction. Intended principally for UI browser.
// The underlying container initially has a capacity of 100.
//
// If much additional functionality is added, should consider containment
// instead of inheritance for std::vector<T>.
//
// Member data:
//
// static G4LogicalVolumeStore* fgInstance
//   - Ptr to the single G4LogicalVolumeStore.

// History:
// 18.04.01 G.Cosmo Migrated to STL vector
// 10.07.95 P.Kent  Initial version

#ifndef G4VLOGICALVOLUMESTORE_HH
#define G4VLOGICALVOLUMESTORE_HH

#include "g4std/vector"

#include "G4LogicalVolume.hh"

class G4LogicalVolumeStore : public G4std::vector<G4LogicalVolume*>
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

//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// G4LogicalVolumeStore
//
// Class description:
//
// Container for all logical volumes, with functionality derived from
// std::vector<T>. The class is a 'singleton', in that only
// one can exist, and access is provided via the static function
// G4LogicalVolumeStore::GetInstance()
//
// All logical volumes should be registered with G4LogicalVolumeStore,
// and removed on their destruction.
// The underlying container initially has a capacity of 100.
// A map indexed by volume names is also recorded for fast search;
// pointers to volumes with same name are stored in buckets.
//
// If much additional functionality is added, should consider containment
// instead of inheritance for std::vector<T>.

// 10.07.95, P.Kent, G.Cosmo - Initial version
// --------------------------------------------------------------------
#ifndef G4LOGICALVOLUMESTORE_HH
#define G4LOGICALVOLUMESTORE_HH 1

#include <vector>
#include <map>

#include "G4LogicalVolume.hh"
#include "G4VStoreNotifier.hh"

class G4LogicalVolumeStore : public std::vector<G4LogicalVolume*>
{
  public:

    static void Register(G4LogicalVolume* pVolume);
      // Add the logical volume to the collection.
    static void DeRegister(G4LogicalVolume* pVolume);
      // Remove the logical volume from the collection.
    static G4LogicalVolumeStore* GetInstance();
      // Get a ptr to the unique G4LogicalVolumeStore, creating it if necessary.
    static void SetNotifier(G4VStoreNotifier* pNotifier);
      // Assign a notifier for allocation/deallocation of the logical volumes.
    static void Clean();
      // Delete all volumes from the store.

    G4LogicalVolume* GetVolume(const G4String& name, G4bool verbose=true,
                  G4bool reverseSearch=false) const;
      // Return the pointer of the first or last volume in the collection having
      // that name. Uses the internal map for fast search and warns if
      // a volume in the collection is not unique or not found.

    inline G4bool IsMapValid() const  { return mvalid; }
    inline void SetMapValid(G4bool val)  { mvalid = val; }
      // Accessor to assess validity of the internal map.
    inline const std::map<G4String,
            std::vector<G4LogicalVolume*> >& GetMap() const { return bmap; }
      // Return the internal map.
    void UpdateMap();
      // Bring contents of internal map up to date and resets validity flag.

    virtual ~G4LogicalVolumeStore();
      // Destructor: takes care to delete allocated logical volumes.

    G4LogicalVolumeStore(const G4LogicalVolumeStore&) = delete;
    G4LogicalVolumeStore& operator=(const G4LogicalVolumeStore&) = delete;
      // Forbidden copy constructor and assignment operator

  protected:

    G4LogicalVolumeStore();

  private:

    static G4LogicalVolumeStore* fgInstance;
    static G4ThreadLocal G4VStoreNotifier* fgNotifier;
    static G4ThreadLocal G4bool locked;

    std::map<G4String, std::vector<G4LogicalVolume*> > bmap;
    G4bool mvalid = false;  // Flag to indicate if map is up to date or not
};

#endif

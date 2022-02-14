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
// G4RegionStore
//
// Class description:
//
// Container for all regions, with functionality derived from
// std::vector<T>. The class is a 'singleton', in that only
// one can exist, and access is provided via the static method
// G4RegionStore::GetInstance().
//
// All regions should be registered with G4RegionStore, and removed on their
// destruction. The underlying container initially has a capacity of 20.
// A map indexed by volume names is also recorded for fast search;
// pointers to regions with same name are stored in buckets.
//
// If much additional functionality is added, should consider containment
// instead of inheritance for std::vector<T>.

// 18.09.02, G.Cosmo - Initial version
// --------------------------------------------------------------------
#ifndef G4REGIONSTORE_HH
#define G4REGIONSTORE_HH 1

#include <vector>
#include <map>

#include "G4Types.hh"
#include "G4String.hh"
#include "G4VStoreNotifier.hh"

class G4Region;
class G4VPhysicalVolume;

class G4RegionStore : public std::vector<G4Region*>
{
  public:

    static void Register(G4Region* pRegion);
      // Add the region to the collection.
    static void DeRegister(G4Region* pRegion);
      // Remove the region from the collection.
    static G4RegionStore* GetInstance();
      // Get a ptr to the unique G4RegionStore, creating it if necessary.
    static void SetNotifier(G4VStoreNotifier* pNotifier);
      // Assign a notifier for allocation/deallocation of regions.
    static void Clean();
      // Delete all regions from the store except for the world region.

    G4bool IsModified() const;
      // Loops through all regions to verify if a region has been
      // modified. It returns TRUE if just one region is modified.
    void ResetRegionModified();
      // Loops through all regions to reset flag for modification
      // to FALSE. Used by the run manager to notify that the
      // physics table has been updated.

    void UpdateMaterialList(G4VPhysicalVolume* currentWorld=nullptr);
      // Forces recomputation of material lists in all regions
      // in the store.

    G4Region* GetRegion(const G4String& name, G4bool verbose = true) const;
      // Returns a region through its name specification. Uses the internal
      // map for fast search and warns if a region in the collection is not
      // unique or not found.

    inline G4bool IsMapValid() const  { return mvalid; }
    inline void SetMapValid(G4bool val)  { mvalid = val; }
      // Accessor to assess validity of the internal map.
    inline const std::map<G4String,
            std::vector<G4Region*> >& GetMap() const { return bmap; }
      // Return the internal map.
    void UpdateMap();
      // Bring contents of internal map up to date and resets validity flag.

    G4Region* FindOrCreateRegion(const G4String& name);
      // Returns a region through its name specification, if it exists.
      // If it does not exist it will allocate one delegating ownership
      // to the client.

    void SetWorldVolume();
      // Set a world volume pointer to a region that belongs to it.
      // Scan over all world volumes.
      // This method should be exclusively used by G4RunManagerKernel.

    virtual ~G4RegionStore();
      // Destructor: takes care to delete allocated regions.

    G4RegionStore(const G4RegionStore&) = delete;
    G4RegionStore& operator=(const G4RegionStore&) = delete;
      // Forbidden copy constructor and assignment operator

  protected:

    G4RegionStore();
      // Protected singleton constructor.

  private:

    static G4RegionStore* fgInstance;
    static G4ThreadLocal G4VStoreNotifier* fgNotifier;
    static G4ThreadLocal G4bool locked;

    std::map<G4String, std::vector<G4Region*> > bmap;
    G4bool mvalid = false;  // Flag to indicate if map is up to date or not
};

#endif

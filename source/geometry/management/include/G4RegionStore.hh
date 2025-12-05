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

// Author: Gabriele Cosmo (CERN), 18.09.2002
// --------------------------------------------------------------------
#ifndef G4REGIONSTORE_HH
#define G4REGIONSTORE_HH

#include <vector>
#include <map>

#include "G4Types.hh"
#include "G4String.hh"
#include "G4VStoreNotifier.hh"

class G4Region;
class G4VPhysicalVolume;

/**
 * @brief G4RegionStore is a singleton class, acting as container
 * for all geometrical regions, with functionality derived from std::vector<T>.
 * All regions should be registered with G4RegionStore, and removed on their
 * destruction. The underlying container initially has a capacity of 20.
 * A map indexed by volume names is also recorded for fast search; pointers
 * to regions with same name are stored in buckets.
*/

class G4RegionStore : public std::vector<G4Region*>
{
  public:

    /**
     * Destructor: takes care to delete allocated regions.
     */
    virtual ~G4RegionStore();

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4RegionStore(const G4RegionStore&) = delete;
    G4RegionStore& operator=(const G4RegionStore&) = delete;

    /**
     * Adds the region 'pRegion' to the collection.
     */
    static void Register(G4Region* pRegion);

    /**
     * Removes the region 'pRegion' from the collection.
     */
    static void DeRegister(G4Region* pRegion);

    /**
     * Returns a pointer to the unique instance of G4RegionStore,
     * creating it if necessary.
     */
    static G4RegionStore* GetInstance();

    /**
     * Assigns a notifier for allocation/deallocation of regions.
     */
    static void SetNotifier(G4VStoreNotifier* pNotifier);

    /**
     * Deletes all regions from the store, except for the world region.
     */
    static void Clean();

    /**
     * Loops through all regions to verify if a region has been modified.
     *  @returns true if just one region is modified.
     */
    G4bool IsModified() const;

    /**
     * Loops through all regions to reset the flag for modification to false.
     * Used by the run manager to notify that the physics table has been updated.
     */
    void ResetRegionModified();

    /**
     * Forces recomputation of the material lists in all regions in the store.
     */
    void UpdateMaterialList(G4VPhysicalVolume* currentWorld = nullptr);

    /**
     * Returns a pointer to a region through its name specification.
     * It uses the internal map for fast search and warns if a region in the
     * collection is not unique or not found.
     */
    G4Region* GetRegion(const G4String& name, G4bool verbose = true) const;

    /**
     * Accessor and modifier to assess validity of the internal map.
     */
    inline G4bool IsMapValid() const { return mvalid; }
    inline void SetMapValid(G4bool val) { mvalid = val; }

    /**
     * Returns the internal map.
     */
    inline const std::map<G4String,
            std::vector<G4Region*> >& GetMap() const { return bmap; }

    /**
     * Brings contents of the internal map up to date and resets validity flag.
     */
    void UpdateMap();

    /**
     * Returns a pointer to a region through its name specification, if it
     * exists. If it does not exist, it will allocate one, delegating ownership
     * to the client.
     */
    G4Region* FindOrCreateRegion(const G4String& name);

    /**
     * Sets a world volume pointer to a region that belongs to it.
     * Scans over all world volumes. This method should be exclusively
     * used by G4RunManagerKernel.
     */
    void SetWorldVolume();

  protected:

    /**
     * Protected singleton constructor.
     */
    G4RegionStore();

  private:

    static G4RegionStore* fgInstance;
    static G4ThreadLocal G4VStoreNotifier* fgNotifier;
    static G4ThreadLocal G4bool locked;

    std::map<G4String, std::vector<G4Region*> > bmap;
    G4bool mvalid = false;  // Flag to indicate if map is up to date or not
};

#endif

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
// G4SolidStore
//
// Class description:
//
// Container for all solids, with functionality derived from
// std::vector<T>. The class is a 'singleton', in that only
// one can exist, and access is provided via the static method
// G4SolidStore::GetInstance().
//
// All solids should be registered with G4SolidStore, and removed on their
// destruction. The underlying container initially has a capacity of 100.
// A map indexed by solid names is also recorded for fast search;
// pointers to solids with same name are stored in buckets.
//
// If much additional functionality is added, should consider containment
// instead of inheritance for std::vector<T>.

// Authors: Gabriele Cosmo & Paul Kent (CERN), 10.07.1995 - Initial version
// --------------------------------------------------------------------
#ifndef G4VSOLIDSTORE_HH
#define G4VSOLIDSTORE_HH

#include <vector>
#include <map>

#include "G4VSolid.hh"
#include "G4VStoreNotifier.hh"

/**
 * @brief G4LogicalVolumeStore is a singleton class, acting as container
 * for all solids primitives, with functionality derived from std::vector<T>.
 * All solids should be registered with G4SolidStore, and removed on their
 * destruction. The underlying container initially has a capacity of 100.
 * A map indexed by solid names is also recorded for fast search; pointers
 * to solids with same name are stored in buckets.
*/

class G4SolidStore : public std::vector<G4VSolid*>
{
  public:

    /**
     * Destructor: takes care to delete allocated solids.
     */
    virtual ~G4SolidStore();

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4SolidStore(const G4SolidStore&) = delete;
    G4SolidStore& operator=(const G4SolidStore&) = delete;

    /**
     * Adds the solid 'pSolid' to the collection.
     */
    static void Register(G4VSolid* pSolid);

    /**
     * Removes the solid 'pSolid' from the collection.
     */
    static void DeRegister(G4VSolid* pSolid);

    /**
     * Returns a pointer to the unique instance of G4SolidStore,
     * creating it if necessary.
     */
    static G4SolidStore* GetInstance();

    /**
     * Assigns a notifier for allocation/deallocation of the solids.
     */
    static void SetNotifier(G4VStoreNotifier* pNotifier);

    /**
     * Deletes all solids from the store.
     */
    static void Clean();

    /**
     * Returns a pointer to the first or last solid in the collection having
     * that name. Uses the internal map for fast search and warns if a solid
     * in the collection is not unique or not found.
     *  @param[in] name The name of the solid to search.
     *  @param[in] verbose Flag for enabling verbosity (default true).
     *  @param[in] reverseSearch Flag to enable inverse search (default false).
     */
    G4VSolid* GetSolid(const G4String& name, G4bool verbose = true,
                       G4bool reverseSearch = false) const;

    /**
     * Accessor and modifier to assess validity of the internal map.
     */
    inline G4bool IsMapValid() const  { return mvalid; }
    inline void SetMapValid(G4bool val)  { mvalid = val; }

    /**
     * Returns the internal map.
     */
    inline const std::map<G4String,
            std::vector<G4VSolid*> >& GetMap() const { return bmap; }

    /**
     * Brings contents of the internal map up to date and resets validity flag.
     */
    void UpdateMap();

  protected:

    /**
     * Protected singleton constructor.
     */
    G4SolidStore();

  private:

    static G4SolidStore* fgInstance;
    static G4ThreadLocal G4VStoreNotifier* fgNotifier;
    static G4ThreadLocal G4bool locked;

    std::map<G4String, std::vector<G4VSolid*> > bmap;
    G4bool mvalid = false;  // Flag to indicate if map is up to date or not
};

#endif

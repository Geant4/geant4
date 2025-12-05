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
// G4Region
//
// Class description:
//
// Defines a region or a group of regions in the detector geometry
// setup, sharing properties associated to materials or production
// cuts which may affect or bias specific physics processes.

// Author: Gabriele Cosmo (CERN), 18.09.2002
// --------------------------------------------------------------------
#ifndef G4REGION_HH
#define G4REGION_HH 1

#include <vector>
#include <map>
#include <algorithm>

#include "G4Types.hh"
#include "G4String.hh"
#include "G4GeomSplitter.hh"

class G4ProductionCuts;
class G4LogicalVolume;
class G4Material;
class G4VUserRegionInformation;
class G4MaterialCutsCouple;
class G4UserLimits;
class G4FieldManager;
class G4FastSimulationManager;
class G4VPhysicalVolume;
class G4UserSteppingAction;

/**
 * @brief G4RegionData encapsulates the fields associated to the class
 * G4Region that may not be read-only..
 */

class G4RegionData
{
  public:

    void initialize()
    {
      fFastSimulationManager = nullptr;
      fRegionalSteppingAction = nullptr;
    }

    G4FastSimulationManager* fFastSimulationManager;
    G4UserSteppingAction* fRegionalSteppingAction;
};

// Encapsulate the methods used by both the master thread and worker threads
// to allocate memory space for the fields encapsulated by the class
// G4RegionData.
//
using G4RegionManager = G4GeomSplitter<G4RegionData>;

/**
 * @brief G4Region defines a region or a group of regions in the detector
 * geometry setup, sharing properties associated to materials or production
 * cuts which may affect or bias specific physics processes.
 */

class G4Region
{
  public:

    /**
     * Constructor for G4Region.
     *  @param[in] name The name of the region.
     */
    G4Region(const G4String& name);

    /**
     * Destructor.
     */
    ~G4Region();

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4Region(const G4Region&) = delete;
    G4Region& operator=(const G4Region&) = delete;

    /**
     * Equality operator, defined by address only.
     */
    inline G4bool operator==(const G4Region& rg) const;

    /**
     * Adds a root logical volume and sets its daughters flags as regions.
     * It also recompute the materials list for the region.
     * Search in the tree can be turned off, assuming the user guarantees
     * the logical volume is NOT already inserted, in which case significant
     * speedup can be achieved in complex flat geometry setups.
     *  @param[in] lv Pointer to the logical volume to act as root region.
     *  @param[in] search To enable/disable search in the tree (default true).
     */
    void AddRootLogicalVolume(G4LogicalVolume* lv, G4bool search=true);

    /**
     * Removes a root logical volume and resets its daughters flags as regions.
     * It also recompute the materials list for the region.
     * The flag for scanning the subtree is always enabled by default.
     *  @param[in] lv Pointer to the logical volume to remove as root region.
     *  @param[in] scan To enable/disable scanning the tree (default true).
     */
    void RemoveRootLogicalVolume(G4LogicalVolume* lv, G4bool scan=true);

    /**
     * Setter/getter for the region's name.
     */
    void SetName(const G4String& name);
    inline const G4String& GetName() const;

    /**
     * Accessors to flag identifying if a region has been modified (and still
     * cuts needs to be computed) or not.
     */
    inline void RegionModified(G4bool flag);
    inline G4bool IsModified() const;

    /**
     * Setter/getter for the production cuts values.
     */
    inline void SetProductionCuts(G4ProductionCuts* cut);
    inline G4ProductionCuts* GetProductionCuts() const;

    /**
     * Methods to return iterators to the lists of root logical volumes
     * and materials.
     */
    inline std::vector<G4LogicalVolume*>::iterator
           GetRootLogicalVolumeIterator();
    inline std::vector<G4Material*>::const_iterator
           GetMaterialIterator() const;

    /**
     * Methods to return the number of elements in the lists of materials and
     * root logical volumes.
     */
    inline std::size_t GetNumberOfMaterials() const;
    inline std::size_t GetNumberOfRootVolumes() const;

    /**
     * Clears the material list and recomputes it, looping through each root
     * logical volume in the region.
     */
    void UpdateMaterialList();

    /**
     * Clears the material list.
     */
    void ClearMaterialList();

    /**
     * Scans recursively the 'lv' logical volume tree, retrieves and places
     * all materials in the list if becoming a region.
     */
    void ScanVolumeTree(G4LogicalVolume* lv, G4bool region);

    /**
     * Setter/getter for the user information data.
     */
    inline void SetUserInformation(G4VUserRegionInformation* ui);
    inline G4VUserRegionInformation* GetUserInformation() const;

    /**
     * Setter/getter for the user-limits associated to a region.
     * Once user-limits are set, it will propagate to the daughter volumes.
     */
    inline void SetUserLimits(G4UserLimits* ul);
    inline G4UserLimits* GetUserLimits() const;

    /**
     * Resets the material-cuts-couples map.
     */
    inline void ClearMap();

    /**
     * Method invoked by G4ProductionCutsTable to register the material-cuts
     * couple pair.
     */
    inline void RegisterMaterialCouplePair(G4Material* mat,
                                           G4MaterialCutsCouple* couple);

    /**
     * Finds a G4MaterialCutsCouple which corresponds to the material 'mat'
     * in the region.
     */
    inline G4MaterialCutsCouple* FindCouple(G4Material* mat);

    /**
     * Setter/getter for the fast-simulation manager. The root logical volume
     * that has the region with G4FastSimulationManager becomes an envelope of
     * fast simulation.
     */
    void SetFastSimulationManager(G4FastSimulationManager* fsm);
    G4FastSimulationManager* GetFastSimulationManager() const;
    
    /**
     * Sets the G4FastSimulationManager pointer to the one for the parent
     * region if it exists, otherwise it sets it to null.
     */
    void ClearFastSimulationManager();

    /**
     * Setter/getter for the field manager. The region with assigned
     * field-manager sets the field to the geometrical area associated with
     * it; priority is anyhow given to local fields eventually set to logical
     * volumes.
     */
    inline void SetFieldManager(G4FieldManager* fm);
    inline G4FieldManager* GetFieldManager() const;

    /**
     * Get method for the world physical volume which the region belongs to.
     * A valid pointer will be assigned by G4RunManagerKernel through the
     * G4RegionStore when the geometry is to be closed. Thus, this pointer
     * may be incorrect at the PreInit and Idle state. If the pointer is null
     * at the proper state, this particular region does not belong to any
     * world (maybe not assigned to any volume, etc.).
     */
    inline G4VPhysicalVolume* GetWorldPhysical() const;

    /**
     * Sets the world physical volume if the region belongs to this world.
     * If 'wp' pointer is null, resets the pointer.
     */
    void SetWorld(G4VPhysicalVolume* wp);

    /**
     * Returns whether the region belongs to the given physical volume 'pv'
     * (recursively scanned to the bottom of the hierarchy).
     */
    G4bool BelongsTo(G4VPhysicalVolume* pv) const;

    /**
     * Returns a region that contains this region. Otherwise null returned.
     *  @param[out] unique Returns true if there is only one parent region
     *              containing the current region.
     *  @returns A pointer to the mother region.
     */
    G4Region* GetParentRegion(G4bool& unique) const;

    /**
     * Setter/getter methods for the regional user stepping action.
     */
    void SetRegionalSteppingAction(G4UserSteppingAction* rusa);
    G4UserSteppingAction* GetRegionalSteppingAction() const;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4Region(__void__&);

    /**
     * Returns the instance ID for multi-threading.
     */
    inline G4int GetInstanceID() const;

    /**
     * Returns the private data instance manager for multi-threading.
     */
    static const G4RegionManager& GetSubInstanceManager();

    /**
     * Clears the memory allocated by the MT sub-instance manager.
     */
    static void Clean();

    /**
     * Utility methods to identify if the region is part of the main mass
     * geometry for tracking or part of a parallel geometry.
     */
    inline void UsedInMassGeometry(G4bool val = true);
    inline void UsedInParallelGeometry(G4bool val = true);
    inline G4bool IsInMassGeometry() const;
    inline G4bool IsInParallelGeometry() const;

  private:

    /**
     * Searchs the specified material 'aMaterial' in the material table and
     * if not present adds it.
     */
    inline void AddMaterial (G4Material* aMaterial);

  private:

    using G4RootLVList = std::vector<G4LogicalVolume*>;
    using G4MaterialList = std::vector<G4Material*>;
    using G4MaterialCouplePair = std::pair<G4Material*, G4MaterialCutsCouple*>;
    using G4MaterialCoupleMap = std::map<G4Material*, G4MaterialCutsCouple*>;

    G4String fName;

    G4RootLVList fRootVolumes;
    G4MaterialList fMaterials;
    G4MaterialCoupleMap fMaterialCoupleMap;

    G4bool fRegionMod = true;
    G4ProductionCuts* fCut = nullptr;

    G4VUserRegionInformation* fUserInfo = nullptr;
    G4UserLimits* fUserLimits = nullptr;
    G4FieldManager* fFieldManager = nullptr;

    G4VPhysicalVolume* fWorldPhys = nullptr;

    G4bool fInMassGeometry = false;
    G4bool fInParallelGeometry = false;

    /** This field is used as instance ID. */
    G4int instanceID;

    /** This field helps to use the class G4RegionManager introduced above. */
    G4GEOM_DLL static G4RegionManager subInstanceManager;
};

#include "G4Region.icc"

// NOTE:
//
// The type G4RegionManager is introduced to encapsulate the methods used by
// both the master thread and worker threads to allocate memory space for
// the fields encapsulated by the class G4RegionData. When each thread
// initializes the value for these fields, it refers to them using a macro
// definition defined below. For every G4Region instance, there is a
// corresponding G4RegionData instance. All G4RegionData instances are
// organized by the class G4RegionManager as an array.
// The field "int instanceID" is added to the class G4Region.
// The value of this field in each G4Region instance is the subscript
// of the corresponding G4RegionData instance.
// In order to use the class G4RegionManager, we add a static member in
// the class G4Region as follows: "static G4RegionManager subInstanceManager".
// For the master thread, the array for G4RegionData instances grows
// dynamically along with G4Region instances are created. For each worker
// thread, it copies the array of G4RegionData instances from the master thread.
// In addition, it invokes a method similiar to the constructor explicitly
// to achieve the partial effect for each instance in the array.

#endif

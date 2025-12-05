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
// G4LogicalVolume
//
// Class description:
//
// Represents a leaf node or unpositioned subtree in the geometry hierarchy.
// Logical volumes are named, and may have daughters ascribed to them.
// They are responsible for retrieval of the physical and tracking attributes
// of the physical volume that it represents: solid, material, magnetic field,
// and optionally, user limits, sensitive detectors, regions, biasing weights.
//
// Get and Set functionality is provided for all attributes, but note that
// most set functions should not be used  when the geometry is `closed'.
// As a  further development, `Guard' checks can be added to ensure
// only legal operations at tracking time.
//
// On construction, solid, material and name must be specified.
//
// Daughters are ascribed and managed by means of a simple
// GetNoDaughters,Get/SetDaughter(n),AddDaughter interface.
//
// Smart voxels as used for tracking optimisation. They're also an attribute.
//
// Logical volumes self register to the logical volume Store on construction,
// and deregister on destruction.
//
// NOTE: This class is NOT meant to act as base class, except for exceptional
//       circumstances of extended types used in the kernel.
//
// Data members:
//
//    std::vector<G4VPhysicalVolume*> fDaughters
//    - Vector of daughters. Given initial size of 0.
//    G4FieldManager* fFieldManager
//    - Pointer (possibly 0) to (magnetic or other) field manager object.
//    G4Material* fMaterial
//    - Pointer to material at this node.
//    G4String fName
//    - Name of logical volume.
//    G4VSensitiveDetector *fSensitiveDetector
//    - Pointer (possibly 0) to `Hit' object.
//    G4VSolid* fSolid
//    - Pointer to solid.
//    G4UserLimits* fUserLimits
//    - Pointer (possibly 0) to user Step limit object for this node.
//    G4SmartVoxelHeader* fVoxel
//    - Pointer (possibly 0) to optimisation info objects.
//    G4bool fOptimise
//    - Flag to identify if optimisation should be applied or not.
//    G4bool fRootRegion
//    - Flag to identify if the logical volume is a root region.
//    G4double fSmartless
//    - Quality for optimisation, average number of voxels to be spent
//      per content.
//    const G4VisAttributes* fVisAttributes
//    - Pointer (possibly 0) to visualization attributes.
//    G4Region* fRegion
//    - Pointer to the cuts region (if any)
//    G4MaterialCutsCouple* fCutsCouple
//    - Pointer (possibly 0) to associated production cuts.
//    G4double fBiasWeight
//    - Weight used in the event biasing technique.
//
// Following data members has been moved to G4Region - M.Asai (Aug/18/2005)
//    G4FastSimulationManager* fFastSimulationManager
//    - Pointer (possibly 0) to G4FastSimulationManager object.
//    G4bool fIsEnvelope
//    - Flags if the Logical Volume is an envelope for a FastSimulationManager.

// Author: Paul Kent (CERN), 11.07.1995 - Initial version
// ------------------------------------------------------------------------
#ifndef G4LOGICALVOLUME_HH
#define G4LOGICALVOLUME_HH

#include <vector>
#include <memory>

#include "G4Types.hh"
#include "G4Region.hh"           // Required by inline methods
#include "G4VPhysicalVolume.hh"  // Need operator == for vector fdaughters
#include "G4GeomSplitter.hh"     // Needed for MT RW data splitting
#include "G4Threading.hh"

// Forward declarations
//
class G4FieldManager;
class G4Material;
class G4VSensitiveDetector;
class G4VSolid;
class G4UserLimits;
class G4SmartVoxelHeader;
class G4FastSimulationManager;
class G4MaterialCutsCouple;
class G4VisAttributes;

/**
 * @brief G4LVData encapsulates the fields associated to the class
 * G4LogicalVolume that may not be read-only. 
 */

class G4LVData
{
  public:

    G4LVData();
    void initialize()
    {
      fSolid = nullptr;
      fSensitiveDetector = nullptr;
      fFieldManager = nullptr;
      fMaterial = nullptr;
      fMass = 0.0;
      fCutsCouple = nullptr;
    }

  public:

    G4VSolid* fSolid = nullptr;
      // Pointer to solid.
    G4VSensitiveDetector* fSensitiveDetector = nullptr;
      // Pointer to sensitive detector.
    G4FieldManager* fFieldManager = nullptr;
      // Pointer (possibly nullptr) to (magnetic or other) field manager object.
    G4Material* fMaterial = nullptr;
      // Pointer to material at this node.
    G4double fMass = 0.0;
      // Mass of the logical volume tree.
    G4MaterialCutsCouple* fCutsCouple = nullptr;
      // Pointer (possibly nullptr) to associated production cuts.
};

/** G4LVManager encapsulates the methods used by both the master thread and
    worker threads to allocate memory space for the fields encapsulated by the
    class G4LVData. */
using G4LVManager = G4GeomSplitter<G4LVData>;

/**
 * @brief G4LogicalVolume represents a leaf node or unpositioned subtree in the
 * geometry hierarchy. Logical volumes are named, and may have daughters
 * ascribed to them. They are responsible for retrieval of the physical and
 * tracking attributes of the physical volume that it represents: solid,
 * material, magnetic field, and optionally, user limits, sensitive detectors,
 * regions, biasing weights.
 */

class G4LogicalVolume
{
  public:
    
    /**
     * Constructor for G4LogicalVolume. The solid and material pointer must be
     * non null. The parameters for field, detector and user limits are optional.
     * The volume also enters itself into the logical volume Store.
     * Optimisation of the geometry (voxelisation) for the volume hierarchy is
     * applied by default. For parameterised volumes in the hierarchy,
     * optimisation is -always- applied.
     *  @param[in] pSolid Pointer to the associated solid primitive.
     *  @param[in] pMaterial Pointer to the associated material.
     *  @param[in] name The volume name.
     *  @param[in] pFieldMgr Pointer to optional magnetic field manager.
     *  @param[in] pSDetector Pointer to optional associated sensitive detector.
     *  @param[in] pULimits Pointer to optional user limits.
     *  @param[in] optimise Flag to enable/disable optimisation structure.
     */
    G4LogicalVolume(G4VSolid* pSolid,
                    G4Material* pMaterial,
              const G4String& name,
                    G4FieldManager* pFieldMgr = nullptr,
                    G4VSensitiveDetector* pSDetector = nullptr,
                    G4UserLimits* pULimits = nullptr,
                    G4bool optimise = true);

    /**
     * Destructor. Removes the logical volume from the logical volume Store.
     * This class is NOT meant to act as base class, except for exceptional
     * circumstances of extended types used in the kernel.
     */
    virtual ~G4LogicalVolume();

    /**
     * Copy-constructor and assignment operator not allowed.
     */
    G4LogicalVolume(const G4LogicalVolume&) = delete;
    G4LogicalVolume& operator=(const G4LogicalVolume&) = delete;

    /**
     * Returns and sets the name of the logical volume.
     */
    inline const G4String& GetName() const;
    void SetName(const G4String& pName);

    /**
     * Returns the number of daughters (0 to n).
     */
    inline std::size_t GetNoDaughters() const;

    /**
     * Returns the ith daughter. Note numbering starts from 0,
     * and no bounds checking is performed.
     */
    inline G4VPhysicalVolume* GetDaughter(const std::size_t i) const;

    /**
     * Adds the volume 'p' as a daughter of the current logical volume.
     */
    void AddDaughter(G4VPhysicalVolume* p);

    /**
     * Returns true if the volume 'p' is a daughter of the current
     * logical volume.
     */
    inline G4bool IsDaughter(const G4VPhysicalVolume* p) const;

    /**
     * Returns true if the volume 'p' is part of the hierarchy of volumes
     * established by the current logical volume. Scans recursively the volume
     * tree.
     */
    G4bool IsAncestor(const G4VPhysicalVolume* p) const;

    /**
     * Removes the volume 'p' from the list of daughters of the current
     * logical volume.
     */
    void RemoveDaughter(const G4VPhysicalVolume* p);

    /**
     * Clears the list of daughters. Used by the physical volume store when
     * the geometry tree is cleared, since modified at run-time.
     */
    void ClearDaughters();

    /**
     * Returns the total number of physical volumes (replicated or placed)
     * in the tree represented by the current logical volume.
     */
    G4int TotalVolumeEntities() const;

    /**
     * Characterises the daughters of this logical volume.
     */
    inline EVolume CharacteriseDaughters() const;

    /**
     * Utility method used by CharacteriseDaughters().
     */
    inline EVolume DeduceDaughtersType() const;
   
    /**
     * Gets and sets the current solid.
     */
    G4VSolid* GetSolid() const;
    void SetSolid(G4VSolid* pSolid);

    /**
     * Gets and sets the current material.
     */
    G4Material* GetMaterial() const;
    void SetMaterial(G4Material* pMaterial);

    /**
     * Sets the material and corresponding Material-Cuts-Couple.
     * This method is invoked by G4Navigator while it is navigating through 
     * material parameterisation.
     */
    void UpdateMaterial(G4Material* pMaterial);

    /**
     * Returns the mass of the logical volume tree computed from the
     * estimated geometrical volume of each solid and material associated
     * to the logical volume and (by default) to its daughters.
     *  @note The computation may require a considerable amount of time,
     *        depending from the complexity of the geometry tree.
     *        The returned value is cached and can be used for successive
     *        calls (default), unless recomputation is forced by providing
     *        'true' for the Boolean argument in input. Computation should
     *        be forced if the geometry setup has changed after the previous
     *        call. By setting the 'propagate' Boolean flag to 'false' the 
     *        method returns the mass of the present logical volume only 
     *        (subtracted for the volume occupied by the daughter volumes).
     *        An optional argument to specify a material is also provided.
     *  @param[in] forced Flag to force recomputation or use cached value.
     *  @param[in] propagate Flag to limit or not computation to daughters.
     *  @param[in] parMaterial Optional pointer to a custom material, usually
     *             used in parameterisations.
     */
    G4double GetMass(G4bool forced = false, G4bool propagate = true,
                     G4Material* parMaterial = nullptr);

    /**
     * Resets the cached value of mass. Ensures that cached value of mass is
     * invalidated due to change in  state, e.g. change of the size of the
     * solid, change of the type of solid, or the addition/deletion of a
     * daughter volume.
     */
    void ResetMass(); 
 
    /**
     * Gets current Field Manager pointer.
     */
    G4FieldManager* GetFieldManager() const;

    /**
     * Sets the Field Manager and propagates it.
     *  @param[in] pFieldMgr Pointer to the field manager.
     *  @param[in] forceToAllDaughters Flag to force propagation to all
     *             daughters (true), or only to daughters having null pointer
     *             to the field manager (false).
     */
    void SetFieldManager(G4FieldManager* pFieldMgr, G4bool forceToAllDaughters); 
 
    /**
     * Gets and sets the current sensitive detector (can be a null pointer).
     */
    G4VSensitiveDetector* GetSensitiveDetector() const;
    void SetSensitiveDetector(G4VSensitiveDetector* pSDetector);

    /**
     * Gets and sets the current User Limits.
     */
    inline G4UserLimits* GetUserLimits() const;
    inline void SetUserLimits(G4UserLimits *pULimits);

    /**
     * Gets and sets the current Voxel Header.
     */
    inline G4SmartVoxelHeader* GetVoxelHeader() const;
    inline void SetVoxelHeader(G4SmartVoxelHeader *pVoxel);
    
    /**
     * Gets and sets the user defined optimisation quality associated to
     * the volume.
     */
    inline G4double GetSmartless() const;
    inline void SetSmartless(G4double s);

    /**
     * Replies if the geometry optimisation (voxelisation) is to be applied
     * for this volume hierarchy.
     */
    inline G4bool IsToOptimise() const;

    /**
     * Specifies if to apply or not geometry optimisation to the volume
     * hierarchy. For parameterised volumes in the hierarchy, optimisation is
     * always applied.
     */
    inline void SetOptimisation(G4bool optim);

    /**
     * Replies if the logical volume represents a root region or not.
     */
    inline G4bool IsRootRegion() const;

    /**
     * Sets/unsets the volume as a root region for cuts.
     */
    inline void SetRegionRootFlag(G4bool rreg);

    /**
     * Replies if the logical volume is part of a cuts region or not.
     */
    inline G4bool IsRegion() const;

    /**
     * Sets/unsets the volume as cuts region.
     */
    inline void SetRegion(G4Region* reg);

    /**
     * Returns the region to which the volume belongs, if any.
     */
    inline G4Region* GetRegion() const;

    /**
     * Propagates region pointer to daughters.
     */
    inline void PropagateRegion();

    /**
     * Accessor and modifier for production cuts.
     */
    const G4MaterialCutsCouple* GetMaterialCutsCouple() const;
    void SetMaterialCutsCouple(G4MaterialCutsCouple* cuts);

    /**
     * Equality defined by address only.
     * Returns true if objects are at same address, else false.
     */
    G4bool operator == (const G4LogicalVolume& lv) const;

    /**
     * Accessor and modifiers for visualization attributes.
     * Arguments are converted to shared_ptr.
     */
    const G4VisAttributes* GetVisAttributes () const;
    void SetVisAttributes (const G4VisAttributes* pVA);
    void SetVisAttributes (const G4VisAttributes& VA);

    /**
     * Gets the current FastSimulationManager pointer if existing,
     * otherwise null.
     */
    inline G4FastSimulationManager* GetFastSimulationManager () const;

    /**
     * Sets and gets the bias weight.
     */
    inline void SetBiasWeight (G4double w);
    inline G4double GetBiasWeight() const;

    /**
     * Returns true if it is not a base-class object.
     */
    virtual G4bool IsExtended() const;

    /**
     * Returns current Field Manager for the master thread.
     */
    inline G4FieldManager* GetMasterFieldManager() const;

    /**
     * Returns current Sensitive Detector for the master thread.
     */
    inline G4VSensitiveDetector* GetMasterSensitiveDetector() const;

    /**
     * Returns current Solid for the master thread.
     */
    inline G4VSolid* GetMasterSolid() const;
  
    /**
     * Returns the instance ID.
     */
    inline G4int GetInstanceID() const;

    /**
     * Returns the private data instance manager.
     */
    static const G4LVManager& GetSubInstanceManager();

    /**
     * Clears memory allocated by sub-instance manager.
     */
    static void Clean();

    /**
     * Sets lock identifier for final deletion of entity.
     */
    inline void Lock();

    /**
     * This method is similar to the constructor. It is used by each worker
     * thread to achieve the partial effect as that of the master thread.
     */
    void InitialiseWorker(G4LogicalVolume* ptrMasterObject,
                          G4VSolid* pSolid, G4VSensitiveDetector* pSDetector);

    /**
     * This method is similar to the destructor. It is used by each worker
     * thread to achieve the partial effect as that of the master thread.
     */
    void TerminateWorker(G4LogicalVolume* ptrMasterObject);

    /**
     * Sets the Field Manager only at this level (does not push down hierarchy)
     */
    void AssignFieldManager(G4FieldManager* fldMgr);
  
    /**
     * Optimised methods, passing thread instance of worker data.
     */
    static G4VSolid* GetSolid(G4LVData& instLVdata);
    static void SetSolid(G4LVData& instLVdata, G4VSolid* pSolid);

    /**
     * Changes the type of the daughters volume to be of type 'atype'.
     * Meant for the user adopting an external navigator for the contents
     * of a volume.
     *  @returns Success (true) or failure (false).
     */
    G4bool ChangeDaughtersType(EVolume atype);

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4LogicalVolume(__void__&);

  private:

    using G4PhysicalVolumeList = std::vector<G4VPhysicalVolume *>;

    /** This field helps in the use of the class G4LVManager. */
    G4GEOM_DLL static G4LVManager subInstanceManager;

    /** Vector of daughters. Given initial size of 0. */
    G4PhysicalVolumeList fDaughters;

    /** Name of logical volume. */
    G4String fName;

    /** Pointer (possibly nullptr) to user step limit object for this node. */
    G4UserLimits* fUserLimits = nullptr;

    /** Pointer (possibly nullptr) to optimisation info objects. */
    G4SmartVoxelHeader* fVoxel = nullptr;

    /** Optimisation quality, average number of voxels to be spent per content. */
    G4double fSmartless = 2.0;

    /** Pointer to the cuts region (if any). */
    G4Region* fRegion = nullptr;

    /** Weight used in the event biasing technique. */
    G4double fBiasWeight = 1.0;

    /** Pointer to visualization attributes. */
    std::shared_ptr<const G4VisAttributes> fVisAttributes;

    // Shadow of master pointers.
    // Each worker thread can access this field from the master thread
    // through these pointers.
    //
    G4VSolid* fSolid = nullptr;
    G4VSensitiveDetector* fSensitiveDetector = nullptr;
    G4FieldManager* fFieldManager = nullptr;
    G4LVData* lvdata = nullptr;  // For use of object persistency

    /** This new field is used as instance ID. */
    G4int instanceID;

    /** Are contents of volume placements, replica, parameterised or external? */
    EVolume fDaughtersVolumeType;

    /** Flag to identify if optimisation should be applied or not. */
    G4bool fOptimise = true;

    /** Flag to identify if the logical volume is a root region. */
    G4bool fRootRegion = false;

    /** Flag to identify if entity is locked for final deletion. */
    G4bool fLock = false;
};

#include "G4LogicalVolume.icc"

// NOTE:
//
// The type G4LVManager is introduced to encapsulate the methods used by
// both the master thread and worker threads to allocate memory space for
// the fields encapsulated by the class G4LVData. When each thread
// initializes the value for these fields, it refers to them using a macro
// definition defined below. For every G4LogicalVolume instance, there is
// a corresponding G4LVData instance. All G4LVData instances are organized
// by the class G4LVManager as an array.
// The field "int instanceID" is added to the class G4LogicalVolume.
// The value of this field in each G4LogicalVolume instance is the subscript
// of the corresponding G4LVData instance.
// In order to use the class G4LVManager, we add a static member in the class
// G4LogicalVolume as follows: "static G4LVManager subInstanceManager".
// For the master thread, the array for G4LVData instances grows dynamically
// along with G4LogicalVolume instances are created. For each worker thread,
// it copies the array of G4LVData instances from the master thread.
// In addition, it invokes a method similiar to the constructor explicitly
// to achieve the partial effect for each instance in the array.

#endif

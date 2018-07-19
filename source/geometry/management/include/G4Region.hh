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
//
// $Id: G4Region.hh 103096 2017-03-15 15:21:33Z gcosmo $
//
// class G4Region
//
// Class description:
//
// Defines a region or a group of regions in the detector geometry
// setup, sharing properties associated to materials or production
// cuts which may affect or bias specific physics processes. 

// History:
// 18.09.02 G.Cosmo Initial version
// --------------------------------------------------------------------
#ifndef G4REGION_HH
#define G4REGION_HH

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

class G4RegionData
{
  // Encapsulates the fields associated to the class
  // G4Region that may not be read-only.

  public:

    void initialize()
    {
      fFastSimulationManager = 0;
      fRegionalSteppingAction = 0;
    }

    G4FastSimulationManager* fFastSimulationManager;
    G4UserSteppingAction* fRegionalSteppingAction;
};

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
//
typedef G4GeomSplitter<G4RegionData> G4RegionManager;

class G4Region
{
  typedef std::vector<G4LogicalVolume*> G4RootLVList;
  typedef std::vector<G4Material*> G4MaterialList;
  typedef std::pair<G4Material*,G4MaterialCutsCouple*> G4MaterialCouplePair;
  typedef std::map<G4Material*,G4MaterialCutsCouple*> G4MaterialCoupleMap;

  public:  // with description

    G4Region(const G4String& name);
    virtual ~G4Region();

    inline G4bool operator==(const G4Region& rg) const;
      // Equality defined by address only.

    void AddRootLogicalVolume(G4LogicalVolume* lv);
    void RemoveRootLogicalVolume(G4LogicalVolume* lv, G4bool scan=true);
      // Add/remove root logical volumes and set/reset their
      // daughters flags as regions. They also recompute the
      // materials list for the region.

    inline void SetName(const G4String& name);
    inline const G4String& GetName() const;
      // Set/get region's name.

    inline void RegionModified(G4bool flag);
    inline G4bool IsModified() const;
      // Accessors to flag identifying if a region has been modified
      // (and still cuts needs to be computed) or not.

    inline void SetProductionCuts(G4ProductionCuts* cut);
    inline G4ProductionCuts* GetProductionCuts() const;

    inline std::vector<G4LogicalVolume*>::iterator
           GetRootLogicalVolumeIterator();
    inline std::vector<G4Material*>::const_iterator
           GetMaterialIterator() const;
      // Return iterators to lists of root logical volumes and materials.

    inline size_t GetNumberOfMaterials() const;
    inline size_t GetNumberOfRootVolumes() const;
      // Return the number of elements in the lists of materials and
      // root logical volumes.

    void UpdateMaterialList();
      // Clears material list and recomputes it looping through
      // each root logical volume in the region.

    void ClearMaterialList();
      // Clears the material list.

    void ScanVolumeTree(G4LogicalVolume* lv, G4bool region);
      // Scans recursively the 'lv' logical volume tree, retrieves
      // and places all materials in the list if becoming a region.

    inline void SetUserInformation(G4VUserRegionInformation* ui);
    inline G4VUserRegionInformation* GetUserInformation() const;
      // Set and Get methods for user information.

    inline void SetUserLimits(G4UserLimits* ul);
    inline G4UserLimits* GetUserLimits() const;
      // Set and Get methods for userL-limits associated to a region.
      // Once user-limits are set, it will propagate to daughter volumes.

    inline void ClearMap();
      // Reset G4MaterialCoupleMap

    inline void RegisterMaterialCouplePair(G4Material* mat,
                                           G4MaterialCutsCouple* couple);
      // Method invoked by G4ProductionCutsTable to register the pair.

    inline G4MaterialCutsCouple* FindCouple(G4Material* mat);
      // Find a G4MaterialCutsCouple which corresponds to the material
      // in this region.

    void SetFastSimulationManager(G4FastSimulationManager* fsm);
    G4FastSimulationManager* GetFastSimulationManager() const;
      // Set and Get methods for G4FastSimulationManager.
      // The root logical volume that has the region with G4FastSimulationManager
      // becomes an envelope of fast simulation.
    
    void ClearFastSimulationManager();
      // Set G4FastSimulationManager pointer to the one for the parent region
      // if it exists. Otherwise set to null.

    inline void SetFieldManager(G4FieldManager* fm);
    inline G4FieldManager* GetFieldManager() const;
      // Set and Get methods for G4FieldManager.
      // The region with assigned field-manager sets the field to the
      // geometrical area associated with it; priority is anyhow given
      // to local fields eventually set to logical volumes.

    inline G4VPhysicalVolume* GetWorldPhysical() const;
      // Get method for the world physical volume which this region
      // belongs to. A valid pointer will be assigned by G4RunManagerKernel
      // through G4RegionStore when the geometry is to be closed. Thus, this
      // pointer may be incorrect at PreInit and Idle state. If the pointer
      // is null at the proper state, this particular region does not belong
      // to any world (maybe not assigned to any volume, etc.).

    void SetWorld(G4VPhysicalVolume* wp);
      // Set the world physical volume if this region belongs to this world.
      // If wp is null, reset the pointer.

    G4bool BelongsTo(G4VPhysicalVolume* thePhys) const;
      // Returns whether this region belongs to the given physical volume
      // (recursively scanned to the bottom of the hierarchy).

    G4Region* GetParentRegion(G4bool& unique) const;
      // Returns a region that contains this region. Otherwise null returned.
      // Flag 'unique' is true if there is only one parent region containing
      // the current region.

    void SetRegionalSteppingAction(G4UserSteppingAction* rusa);
    G4UserSteppingAction* GetRegionalSteppingAction() const;
      // Set/Get method of the regional user stepping action

  public:  // without description

    G4Region(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    inline G4int GetInstanceID() const;
      // Returns the instance ID.

    static const G4RegionManager& GetSubInstanceManager();
      // Returns the private data instance manager. 

    static void Clean();
      // Clear memory allocated by sub-instance manager.

    inline void UsedInMassGeometry(G4bool val=true);
    inline void UsedInParallelGeometry(G4bool val=true);
    inline G4bool IsInMassGeometry() const;
    inline G4bool IsInParallelGeometry() const;
      // Utility methods to identify if region is part of the main mass
      // geometry for tracking or a parallel geometry.

  private:

    G4Region(const G4Region&);
    G4Region& operator=(const G4Region&);
      // Private copy constructor and assignment operator.

    inline void AddMaterial (G4Material* aMaterial);
      // Searchs the specified material in the material table and
      // if not present adds it.

  private:

    G4String fName;

    G4RootLVList fRootVolumes;
    G4MaterialList fMaterials;
    G4MaterialCoupleMap fMaterialCoupleMap;

    G4bool fRegionMod;
    G4ProductionCuts* fCut;

    G4VUserRegionInformation* fUserInfo;
    G4UserLimits* fUserLimits;
    G4FieldManager* fFieldManager;

    G4VPhysicalVolume* fWorldPhys;

    G4bool fInMassGeometry;
    G4bool fInParallelGeometry;

    G4int instanceID;
      // This field is used as instance ID.
    G4GEOM_DLL static G4RegionManager subInstanceManager;
      // This field helps to use the class G4RegionManager introduced above.
};

#include "G4Region.icc"

#endif

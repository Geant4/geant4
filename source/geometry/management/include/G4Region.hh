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
// $Id: G4Region.hh,v 1.14 2005/11/09 14:54:03 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
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

class G4ProductionCuts;
class G4LogicalVolume;
class G4Material;
class G4VUserRegionInformation;
class G4MaterialCutsCouple;
class G4UserLimits;
class G4FastSimulationManager;
class G4VPhysicalVolume;

#include <vector>
#include <map>

#include "G4Types.hh"
#include "G4String.hh"

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
    void RemoveRootLogicalVolume(G4LogicalVolume* lv);
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

    inline void SetFastSimulationManager(G4FastSimulationManager* fsm);
    inline G4FastSimulationManager* GetFastSimulationManager() const;
      // Set and Get methods for G4FastSimulationManager.
      // The root logical volume that has the region with G4FastSimulationManager
      // becomes an envelope of fast simulation.
    
    void ClearFastSimulationManager();
      // Set G4FastSimulationManager pointer to the one for the parent region
      // if it exists. Otherwise set to null.

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
      // (recursively scanned to the bottom of the hierarchy)

    G4Region* GetParentRegion() const;
      // Returns a region that contains this region. Otherwise null returned.

  public:  // without description

    G4Region(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

  private:

    G4Region(const G4Region&);
    G4Region& operator=(const G4Region&);
      // Private copy constructor and assignment operator.

  private:

    G4String fName;

    G4RootLVList fRootVolumes;
    G4MaterialList fMaterials;
    G4MaterialCoupleMap fMaterialCoupleMap;

    G4bool fRegionMod;
    G4ProductionCuts* fCut;

    G4VUserRegionInformation* fUserInfo;
    G4UserLimits* fUserLimits;

    G4FastSimulationManager* fFastSimulationManager;

    G4VPhysicalVolume* fWorldPhys;
};

#include "G4Region.icc"

#endif

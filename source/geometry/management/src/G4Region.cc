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
// G4Region class implementation
//
// 18.09.02, G.Cosmo - Initial version
// --------------------------------------------------------------------

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VNestedParameterisation.hh"
#include "G4VUserRegionInformation.hh"
#include "G4Material.hh"

// These macros changes the references to fields that are now encapsulated
// in the class G4RegionData.
//
#define G4MT_fsmanager ((subInstanceManager.offset[instanceID]).fFastSimulationManager)
#define G4MT_rsaction ((subInstanceManager.offset[instanceID]).fRegionalSteppingAction)

// This new field helps to use the class G4RegionManager
//
G4RegionManager G4Region::subInstanceManager;

// *******************************************************************
// GetSubInstanceManager:
//  - Returns the private data instance manager.
// *******************************************************************
//
const G4RegionManager& G4Region::GetSubInstanceManager()
{
  return subInstanceManager;
}

// *******************************************************************
// Constructor:
//  - Adds self to region Store
// *******************************************************************
//
G4Region::G4Region(const G4String& pName)
  : fName(pName)
{

  instanceID = subInstanceManager.CreateSubInstance();
  G4MT_fsmanager = nullptr;
  G4MT_rsaction = nullptr;

  G4RegionStore* rStore = G4RegionStore::GetInstance();
  if (rStore->GetRegion(pName, false))
  {
    std::ostringstream message;
    message << "The region has NOT been registered !" << G4endl
            << "          Region " << pName << " already existing in store !"
            << G4endl;
    G4Exception("G4Region::G4Region()", "GeomMgt1001",
                JustWarning, message);
  }
  else
  {
    rStore->Register(this);
  }
}

// ********************************************************************
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
// ********************************************************************
//
G4Region::G4Region( __void__& )
  : fName("")
{
  instanceID = subInstanceManager.CreateSubInstance();
  G4MT_fsmanager = nullptr;
  G4MT_rsaction = nullptr;

  // Register to store
  //
  G4RegionStore::GetInstance()->Register(this);
}

// *******************************************************************
// Destructor:
//  - Removes self from region Store
// *******************************************************************
//
G4Region::~G4Region()
{
  G4RegionStore::GetInstance()->DeRegister(this);
  if(fUserInfo != nullptr)  { delete fUserInfo; }
}

// ********************************************************************
// SetName - Set region name and notify store of the change
// ********************************************************************
//
void G4Region::SetName(const G4String& pName)
{
  fName = pName;
  G4RegionStore::GetInstance()->SetMapValid(false);
}

// ********************************************************************
// SetFastSimulationManager
// ********************************************************************
//
void G4Region::SetFastSimulationManager(G4FastSimulationManager* fsm)
{
  G4MT_fsmanager = fsm;
}

// ********************************************************************
// GetFastSimulationManager
// ********************************************************************
//
G4FastSimulationManager* G4Region::GetFastSimulationManager() const
{
  return G4MT_fsmanager;
}

// ********************************************************************
// SetRegionalSteppingAction
// ********************************************************************
//
void G4Region::SetRegionalSteppingAction(G4UserSteppingAction* rusa)
{
  G4MT_rsaction = rusa;
}

// ********************************************************************
// GetRegionalSteppingAction
// ********************************************************************
//
G4UserSteppingAction* G4Region::GetRegionalSteppingAction() const
{
  return G4MT_rsaction;
}

// *******************************************************************
// ScanVolumeTree:
//  - Scans recursively the 'lv' logical volume tree, retrieves
//    and places all materials in the list.
//  - The boolean flag 'region' identifies if the volume tree must
//    have region reset (false) or if the current region must be
//    associated to the logical volume 'lv' and its tree (true).
// *******************************************************************
//
void G4Region::ScanVolumeTree(G4LogicalVolume* lv, G4bool region)
{
  // If logical volume is going to become a region, add 
  // its material to the list if not already present
  //
  G4Region* currentRegion = nullptr;
  std::size_t noDaughters = lv->GetNoDaughters();
  G4Material* volMat = lv->GetMaterial();
  if((volMat == nullptr) && fInMassGeometry)
  {
    std::ostringstream message;
    message << "Logical volume <" << lv->GetName() << ">" << G4endl
            << "does not have a valid material pointer." << G4endl
            << "A logical volume belonging to the (tracking) world volume "
            << "must have a valid material.";
    G4Exception("G4Region::ScanVolumeTree()", "GeomMgt0002",
                FatalException, message, "Check your geometry construction.");
  }
  if (region)
  {
    currentRegion = this;
    if (volMat != nullptr)
    { 
      AddMaterial(volMat); 
      G4Material* baseMat = const_cast<G4Material*>(volMat->GetBaseMaterial());
      if (baseMat != nullptr) { AddMaterial(baseMat); }
    }
  }

  // Set the LV region to be either the current region or NULL,
  // according to the boolean selector
  //
  lv->SetRegion(currentRegion);

  // Stop recursion here if no further daughters are involved
  //
  if(noDaughters==0) return;

  G4VPhysicalVolume* daughterPVol = lv->GetDaughter(0);
  if (daughterPVol->IsParameterised())
  {
    // Adopt special treatment in case of parameterised volumes,
    // where parameterisation involves a new material scan
    //
    G4VPVParameterisation* pParam = daughterPVol->GetParameterisation();

    if (pParam->GetMaterialScanner() != nullptr)
    {
      std::size_t matNo = pParam->GetMaterialScanner()->GetNumberOfMaterials();
      for (std::size_t mat=0; mat<matNo; ++mat)
      {
        volMat = pParam->GetMaterialScanner()->GetMaterial((G4int)mat);
        if(!volMat && fInMassGeometry)
        {
          std::ostringstream message;
          message << "The parameterisation for the physical volume <"
                  << daughterPVol->GetName() << ">" << G4endl
                  << "does not return a valid material pointer." << G4endl
                  << "A volume belonging to the (tracking) world volume must "
                  << "have a valid material.";
          G4Exception("G4Region::ScanVolumeTree()", "GeomMgt0002",
                      FatalException, message, "Check your parameterisation.");
        }
        if (volMat != nullptr)
        { 
          AddMaterial(volMat); 
          G4Material* baseMat = const_cast<G4Material*>(volMat->GetBaseMaterial());
          if (baseMat != nullptr) { AddMaterial(baseMat); }
        }
      }
    }
    else
    {
      std::size_t repNo = daughterPVol->GetMultiplicity();
      for (std::size_t rep=0; rep<repNo; ++rep)
      {
        volMat = pParam->ComputeMaterial((G4int)rep, daughterPVol);
        if((volMat == nullptr) && fInMassGeometry)
        {
          std::ostringstream message;
          message << "The parameterisation for the physical volume <"
                  << daughterPVol->GetName() << ">" << G4endl
                  << "does not return a valid material pointer." << G4endl
                  << "A volume belonging to the (tracking) world volume must "
                  << "have a valid material.";
          G4Exception("G4Region::ScanVolumeTree()", "GeomMgt0002",
                      FatalException, message, "Check your parameterisation.");
        }
        if(volMat != nullptr)
        { 
          AddMaterial(volMat);
          G4Material* baseMat = const_cast<G4Material*>(volMat->GetBaseMaterial());
          if (baseMat != nullptr) { AddMaterial(baseMat); }
        }
      }
    }
    G4LogicalVolume* daughterLVol = daughterPVol->GetLogicalVolume();
    ScanVolumeTree(daughterLVol, region);
  }
  else
  {
    for (std::size_t i=0; i<noDaughters; ++i)
    {
      G4LogicalVolume* daughterLVol = lv->GetDaughter(i)->GetLogicalVolume();
      if (!daughterLVol->IsRootRegion())
      {
        // Set daughter's LV to be a region and store materials in
        // the materials list, if the LV is not already a root region
        //
        ScanVolumeTree(daughterLVol, region);
      }
    }
  }
}

// *******************************************************************
// AddRootLogicalVolume:
//  - Adds a root logical volume and sets its daughters flags as
//    regions. It also recomputes the materials list for the region.
// *******************************************************************
//
void G4Region::AddRootLogicalVolume(G4LogicalVolume* lv, G4bool search)
{
  // Check the logical volume is not already in the list
  //
  if (search) 
  {
    auto pos = std::find(fRootVolumes.cbegin(),fRootVolumes.cend(),lv);
    if (pos == fRootVolumes.cend())
    {
      // Insert the root volume in the list and set it as root region
      //
      fRootVolumes.push_back(lv);
      lv->SetRegionRootFlag(true);
    }
  }
  else  // WARNING: user *MUST* guarantee lv is not already inserted.
  {     // Providing speedup for very complex flat geometries
    fRootVolumes.push_back(lv);
    lv->SetRegionRootFlag(true);
  }
  // Scan recursively the tree of daugther volumes and set regions
  //
  ScanVolumeTree(lv, true);

  // Set region as modified
  //
  fRegionMod = true;
}

// *******************************************************************
// RemoveRootLogicalVolume:
//  - Removes a root logical volume and resets its daughters flags as
//    regions. It also recomputes the materials list for the region.
// *******************************************************************
//
void G4Region::RemoveRootLogicalVolume(G4LogicalVolume* lv, G4bool scan)
{
  // Find and remove logical volume from the list
  //
  auto pos = std::find(fRootVolumes.cbegin(),fRootVolumes.cend(),lv);
  if (pos != fRootVolumes.cend())
  {
    if (fRootVolumes.size() != 1)  // Avoid resetting flag for world since
    {                              // volume may be already deleted !
      lv->SetRegionRootFlag(false);
    }
    fRootVolumes.erase(pos);
  }

  if (scan)  // Update the materials list
  {
    UpdateMaterialList();
  }

  // Set region as modified
  //
  fRegionMod = true;
}

// ********************************************************************
// Clean
// ********************************************************************
//
void G4Region::Clean()
{
  subInstanceManager.FreeSlave();
}

// *******************************************************************
// ClearMaterialList:
//  - Clears the material list.
// *******************************************************************
//
void G4Region::ClearMaterialList()
{
  fMaterials.clear();
}

// *******************************************************************
// UpdateMaterialList:
//  - computes material list looping through
//    each root logical volume in the region.
// *******************************************************************
//
void G4Region::UpdateMaterialList()
{
  // Reset the materials list
  //
  ClearMaterialList();

  // Loop over the root logical volumes and rebuild the list
  // of materials from scratch
  //
  for (auto pLV=fRootVolumes.cbegin(); pLV!=fRootVolumes.cend(); ++pLV)
  {
    ScanVolumeTree(*pLV, true);
  }
}

// *******************************************************************
// SetWorld:
//  - Set the world physical volume if this region belongs to this
//    world. If the given pointer is null, reset the pointer.
// *******************************************************************
//
void G4Region::SetWorld(G4VPhysicalVolume* wp)
{
  if(!wp)
  { fWorldPhys = nullptr; }
  else
  { if(BelongsTo(wp)) fWorldPhys = wp; }

  return;
}

// *******************************************************************
// BelongsTo:
//  - Returns whether this region belongs to the given physical volume
//    (recursively scanned to the bottom of the hierarchy)
// *******************************************************************
// 
G4bool G4Region::BelongsTo(G4VPhysicalVolume* thePhys) const
{
  G4LogicalVolume* currLog = thePhys->GetLogicalVolume();
  if (currLog->GetRegion()==this) {return true;}

  std::size_t nDaughters = currLog->GetNoDaughters();
  while (nDaughters--)  // Loop checking, 06.08.2015, G.Cosmo
  {
    if (BelongsTo(currLog->GetDaughter(nDaughters))) {return true;}
  }

  return false;
}

// *******************************************************************
// ClearFastSimulationManager:
//  - Set G4FastSimulationManager pointer to the one for the parent region
//    if it exists. Otherwise set to null.
// *******************************************************************
//
void G4Region::ClearFastSimulationManager()
{
  G4bool isUnique;
  G4Region* parent = GetParentRegion(isUnique);
  if(parent != nullptr)
  {
    if (isUnique)
    {
      G4MT_fsmanager = parent->GetFastSimulationManager();
    }
    else
    {
      std::ostringstream message;
      message << "Region <" << fName << "> belongs to more than"
              << " one parent region !" << G4endl
              << "A region cannot belong to more than one direct parent region,"
              << G4endl
              << "to have fast-simulation assigned.";
      G4Exception("G4Region::ClearFastSimulationManager()",
                  "GeomMgt1002", JustWarning, message);
      G4MT_fsmanager = nullptr;
    }
  }
  else
  {
    G4MT_fsmanager = nullptr;
  }
}

// *******************************************************************
// GetParentRegion:
//  - Returns a region that contains this region.
//    Otherwise null is returned.
// *******************************************************************
// 
G4Region* G4Region::GetParentRegion(G4bool& unique) const
{
  G4Region* parent = nullptr; unique = true;
  G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();

  // Loop over all logical volumes in the store
  //
  for(auto lvItr=lvStore->cbegin(); lvItr!=lvStore->cend(); ++lvItr)
  {
    std::size_t nD = (*lvItr)->GetNoDaughters();
    G4Region* aR = (*lvItr)->GetRegion();

    // Loop over all daughters of each logical volume
    //
    for(std::size_t iD=0; iD<nD; ++iD)
    {
      if((*lvItr)->GetDaughter(iD)->GetLogicalVolume()->GetRegion()==this)
      { 
        if(parent)
        {
          if(parent!=aR) { unique = false; }
        }
        else  // Cache LV parent region which includes a daughter volume
              // with the same associated region as the current one
        {
          parent = aR;
        }
      }
    }
  }
  return parent;
}

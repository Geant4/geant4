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
// class G4LogicalVolume implementation
//
// 15.01.13 G.Cosmo, A.Dotti: Modified for thread-safety for MT
// 01.03.05 G.Santin: Added flag for optional propagation of GetMass()
// 17.05.02 G.Cosmo: Added flag for optional optimisation
// 12.02.99 S.Giani: Default initialization of voxelization quality
// 04.08.97 P.M.DeFreitas: Added methods for parameterised simulation 
// 11.07.95 P.Kent: Initial version
// --------------------------------------------------------------------

#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VSolid.hh"
#include "G4Material.hh"
#include "G4VPVParameterisation.hh"
#include "G4VisAttributes.hh"

#include "G4UnitsTable.hh"

G4LVData::G4LVData() {;}

// This new field helps to use the class G4LVManager
//
G4LVManager G4LogicalVolume::subInstanceManager;

// These macros change the references to fields that are now encapsulated
// in the class G4LVData.
//
#define G4MT_solid     ((subInstanceManager.offset[instanceID]).fSolid)
#define G4MT_sdetector ((subInstanceManager.offset[instanceID]).fSensitiveDetector)
#define G4MT_fmanager  ((subInstanceManager.offset[instanceID]).fFieldManager)
#define G4MT_material  ((subInstanceManager.offset[instanceID]).fMaterial)
#define G4MT_mass      ((subInstanceManager.offset[instanceID]).fMass)
#define G4MT_ccouple   ((subInstanceManager.offset[instanceID]).fCutsCouple)
#define G4MT_instance  (subInstanceManager.offset[instanceID])

// ********************************************************************
// Constructor - sets member data and adds to logical Store,
//               voxel pointer for optimisation set to 0 by default.
//               Initialises daughter vector to 0 length.
// ********************************************************************
//
G4LogicalVolume::G4LogicalVolume( G4VSolid* pSolid,
                                  G4Material* pMaterial,
                            const G4String& name,
                                  G4FieldManager* pFieldMgr,
                                  G4VSensitiveDetector* pSDetector,
                                  G4UserLimits* pULimits,
                                  G4bool optimise )
 : fDaughters(0,(G4VPhysicalVolume*)nullptr), fDaughtersVolumeType(kNormal),
   fOptimise(optimise)
{
  // Initialize 'Shadow'/master pointers - for use in copying to workers
  //
  fSolid = pSolid;
  fSensitiveDetector = pSDetector;
  fFieldManager = pFieldMgr;

  instanceID = subInstanceManager.CreateSubInstance();
  AssignFieldManager(pFieldMgr);
  
  G4MT_mass = 0.;
  G4MT_ccouple = nullptr;

  SetSolid(pSolid);
  SetMaterial(pMaterial);
  SetName(name);
  SetSensitiveDetector(pSDetector);
  SetUserLimits(pULimits);    

  // Initialize 'Shadow' data structure - for use by object persistency
  //
  lvdata = new G4LVData();
  lvdata->fSolid = pSolid;
  lvdata->fMaterial = pMaterial;

  //
  // Add to store
  //
  G4LogicalVolumeStore::Register(this);
}

// ********************************************************************
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
// ********************************************************************
//
G4LogicalVolume::G4LogicalVolume( __void__& )
 : fDaughters(0,(G4VPhysicalVolume*)nullptr), fName("")
{
  instanceID = subInstanceManager.CreateSubInstance();
  
  SetSensitiveDetector(nullptr);    // G4MT_sdetector = nullptr;
  SetFieldManager(nullptr, false);  // G4MT_fmanager = nullptr;

  G4MT_mass = 0.;
  G4MT_ccouple = nullptr;
  
  // Add to store
  //
  G4LogicalVolumeStore::Register(this);
}

// ********************************************************************
// Destructor - Removes itself from solid Store
// NOTE: Not virtual
// ********************************************************************
//
G4LogicalVolume::~G4LogicalVolume()
{
  if (!fLock && fRootRegion)  // De-register root region first if not locked
  {                           // and flagged as root logical-volume
    fRegion->RemoveRootLogicalVolume(this, true);
  }
  delete lvdata;
  G4LogicalVolumeStore::DeRegister(this);
}

// ********************************************************************
// SetName - Set volume name and notify store of the change
// ********************************************************************
//
void G4LogicalVolume::SetName(const G4String& pName)
{
  fName = pName;
  G4LogicalVolumeStore::GetInstance()->SetMapValid(false);
}

// ********************************************************************
// InitialiseWorker
//
// This method is similar to the constructor. It is used by each worker
// thread to achieve the same effect as that of the master thread exept
// to register the new created instance. This method is invoked explicitly.
// It does not create a new G4LogicalVolume instance. It only assign the value
// for the fields encapsulated by the class G4LVData.
// ********************************************************************
//
void G4LogicalVolume::
InitialiseWorker( G4LogicalVolume* /*pMasterObject*/,
                  G4VSolid* pSolid,
                  G4VSensitiveDetector* pSDetector)
{
  subInstanceManager.SlaveCopySubInstanceArray();

  SetSolid(pSolid);
  SetSensitiveDetector(pSDetector); //  How this object is available now ?
  AssignFieldManager(fFieldManager);
   // Should be set - but a per-thread copy is not available yet
   // Must not call SetFieldManager(), which propagates FieldMgr

#ifdef CLONE_FIELD_MGR
  // Create a field FieldManager by cloning
  //
  G4FieldManager workerFldMgr = fFieldManager->GetWorkerClone(G4bool* created);
  if( created || (GetFieldManager() != workerFldMgr) )
  {
    SetFieldManager(fFieldManager, false); // which propagates FieldMgr
  }
  else
  {
    // Field manager existed and is equal to current one
    //
    AssignFieldManager(workerFldMgr);
  }
#endif
}

// ********************************************************************
// Clean
// ********************************************************************
//
void G4LogicalVolume::Clean()
{
  subInstanceManager.FreeSlave();
}

// ********************************************************************
// TerminateWorker
//
// This method is similar to the destructor. It is used by each worker
// thread to achieve the partial effect as that of the master thread.
// For G4LogicalVolume instances, nothing more to do here.
// ********************************************************************
//
void G4LogicalVolume::
TerminateWorker( G4LogicalVolume* /*pMasterObject*/)
{
}

// ********************************************************************
// GetSubInstanceManager
//
// Returns the private data instance manager.    
// ********************************************************************
//
const G4LVManager& G4LogicalVolume::GetSubInstanceManager()
{
  return subInstanceManager;
}

// ********************************************************************
// GetFieldManager
// ********************************************************************
//
G4FieldManager* G4LogicalVolume::GetFieldManager() const
{
  return G4MT_fmanager;
}

// ********************************************************************
// AssignFieldManager
// ********************************************************************
//
void G4LogicalVolume::AssignFieldManager( G4FieldManager *fldMgr)
{
  G4MT_fmanager= fldMgr;
  if(G4Threading::IsMasterThread())  { fFieldManager = fldMgr; }
}

// ********************************************************************
// IsExtended
// ********************************************************************
//
G4bool G4LogicalVolume::IsExtended() const
{
  return false;
}

// ********************************************************************
// SetFieldManager
// ********************************************************************
//
void
G4LogicalVolume::SetFieldManager(G4FieldManager* pNewFieldMgr,
                                 G4bool          forceAllDaughters) 
{
  AssignFieldManager(pNewFieldMgr);

  auto NoDaughters = GetNoDaughters();
  while ( (NoDaughters--)>0 )
  {
    G4LogicalVolume* DaughterLogVol; 
    DaughterLogVol = GetDaughter(NoDaughters)->GetLogicalVolume();
    if ( forceAllDaughters || (DaughterLogVol->GetFieldManager() == nullptr) )
    {
      DaughterLogVol->SetFieldManager(pNewFieldMgr, forceAllDaughters);
    }
  }
}

// ********************************************************************
// AddDaughter
// ********************************************************************
//
void G4LogicalVolume::AddDaughter(G4VPhysicalVolume* pNewDaughter)
{
  EVolume daughterType = pNewDaughter->VolumeType();

  // The type of the navigation needed is determined by the first daughter
  //
  if( fDaughters.empty() )
  {
    fDaughtersVolumeType = daughterType;
  }
  else
  {
    // Check consistency of detector description

    // 1. A replica or parameterised volume can have only one daughter
    //
    if( fDaughters[0]->IsReplicated() )
    {
      std::ostringstream message;
      message << "ERROR - Attempt to place a volume in a mother volume"
        << G4endl
        << "        already containing a replicated volume." << G4endl
        << "        A volume can either contain several placements" << G4endl
        << "        or a unique replica or parameterised volume !" << G4endl
        << "           Mother logical volume: " << GetName() << G4endl
        << "           Placing volume: " << pNewDaughter->GetName()
        << G4endl;
      G4Exception("G4LogicalVolume::AddDaughter()", "GeomMgt0002",
                  FatalException, message,
                  "Replica or parameterised volume must be the only daughter!");
    }
    else
    {
      // 2. Ensure that Placement and External physical volumes do not mix
      //
      if(  daughterType != fDaughtersVolumeType )
      {
        std::ostringstream message;
        message << "ERROR - Attempt to place a volume in a mother volume"
          << G4endl
          << "        already containing a different type of volume." << G4endl
          << "        A volume can either contain" << G4endl
          << "        - one or more placements, OR" << G4endl
          << "        - one or more 'external' type physical volumes." << G4endl
          << "          Mother logical volume: " << GetName() << G4endl
          << "          Volume being placed: " << pNewDaughter->GetName()
          << G4endl;
        G4Exception("G4LogicalVolume::AddDaughter()", "GeomMgt0002",
                    FatalException, message,
                    "Cannot mix placements and external physical volumes !");
      }
    }
  }

  // Invalidate previous calculation of mass - if any - for all threads
  //
  G4MT_mass = 0.;
  fDaughters.push_back(pNewDaughter);

  G4LogicalVolume* pDaughterLogical = pNewDaughter->GetLogicalVolume();

  // Propagate the Field Manager, if the daughter has no field Manager
  //
  G4FieldManager* pDaughterFieldManager = pDaughterLogical->GetFieldManager();

  // Avoid propagating the fieldManager pointer if null
  // and daughter's one is null as well...
  // 
  if( (G4MT_fmanager != nullptr ) && (pDaughterFieldManager == nullptr) )
  {
    pDaughterLogical->SetFieldManager(G4MT_fmanager, false);
  }
  if (fRegion)
  {
    PropagateRegion();
    fRegion->RegionModified(true);
  }
}

// ********************************************************************
// RemoveDaughter
// ********************************************************************
//
void G4LogicalVolume::RemoveDaughter(const G4VPhysicalVolume* p)
{
  for (auto i=fDaughters.cbegin(); i!=fDaughters.cend(); ++i )
  {
    if (**i==*p)
    {
      fDaughters.erase(i);
      break;
    }
  }
  if (fRegion)
  {
    fRegion->RegionModified(true);
  }
  G4MT_mass = 0.;
}

// ********************************************************************
// ClearDaughters
// ********************************************************************
//
void G4LogicalVolume::ClearDaughters()
{
  fDaughters.erase(fDaughters.cbegin(), fDaughters.cend());
  if (fRegion != nullptr)
  {
    fRegion->RegionModified(true);
  }
  G4MT_mass = 0.;
}

// ********************************************************************
// ResetMass
// ********************************************************************
//
void G4LogicalVolume::ResetMass()
{
  G4MT_mass= 0.0;
}

// ********************************************************************
// GetSolid
// ********************************************************************
//
G4VSolid* G4LogicalVolume::GetSolid(G4LVData &instLVdata) // const
{
  return instLVdata.fSolid;
}

G4VSolid* G4LogicalVolume::GetSolid() const
{
  return this->GetSolid( subInstanceManager.offset[instanceID] ); 
}

// ********************************************************************
// SetSolid
// ********************************************************************
//
void G4LogicalVolume::SetSolid(G4VSolid *pSolid)
{

  G4MT_solid = pSolid;
  this->ResetMass(); 
}

void G4LogicalVolume::SetSolid(G4LVData& instLVdata, G4VSolid* pSolid)
{  
  instLVdata.fSolid = pSolid;
  instLVdata.fMass = 0.0;
}

// ********************************************************************
// GetMaterial
// ********************************************************************
//
G4Material* G4LogicalVolume::GetMaterial() const
{
  return G4MT_material;
}

// ********************************************************************
// SetMaterial
// ********************************************************************
//
void G4LogicalVolume::SetMaterial(G4Material* pMaterial)
{
  G4MT_material = pMaterial;
  G4MT_mass = 0.0;
}

// ********************************************************************
// UpdateMaterial
// ********************************************************************
//
void G4LogicalVolume::UpdateMaterial(G4Material* pMaterial)
{
  G4MT_material=pMaterial;
  if (fRegion != nullptr) { G4MT_ccouple = fRegion->FindCouple(pMaterial); }
  G4MT_mass = 0.0;
}

// ********************************************************************
// GetSensitiveDetector
// ********************************************************************
//
G4VSensitiveDetector* G4LogicalVolume::GetSensitiveDetector() const
{
  return G4MT_sdetector;
}

// ********************************************************************
// SetSensitiveDetector
// ********************************************************************
//
void G4LogicalVolume::SetSensitiveDetector(G4VSensitiveDetector* pSDetector)
{
  G4MT_sdetector = pSDetector;
  if (G4Threading::IsMasterThread())  { fSensitiveDetector = pSDetector; }
}

// ********************************************************************
// GetMaterialCutsCouple
// ********************************************************************
//
const G4MaterialCutsCouple* G4LogicalVolume::GetMaterialCutsCouple() const
{
  return G4MT_ccouple;
}

// ********************************************************************
// SetMaterialCutsCouple
// ********************************************************************
//
void G4LogicalVolume::SetMaterialCutsCouple(G4MaterialCutsCouple* cuts)
{
  G4MT_ccouple = cuts;
}

// ********************************************************************
// IsAncestor
//
// Finds out if the current logical volume is an ancestor of a given 
// physical volume
// ********************************************************************
//
G4bool
G4LogicalVolume::IsAncestor(const G4VPhysicalVolume* aVolume) const
{
  G4bool isDaughter = IsDaughter(aVolume);
  if (!isDaughter)
  {
    for (auto itDau = fDaughters.cbegin(); itDau != fDaughters.cend(); ++itDau)
    {
      isDaughter = (*itDau)->GetLogicalVolume()->IsAncestor(aVolume);
      if (isDaughter)  break;
    }
  }
  return isDaughter;
}

// ********************************************************************
// TotalVolumeEntities
//
// Returns the total number of physical volumes (replicated or placed)
// in the tree represented by the current logical volume.
// ********************************************************************
//
G4int G4LogicalVolume::TotalVolumeEntities() const
{
  G4int vols = 1;
  for (auto itDau = fDaughters.cbegin(); itDau != fDaughters.cend(); ++itDau)
  {
    G4VPhysicalVolume* physDaughter = (*itDau);
    vols += physDaughter->GetMultiplicity()
           *physDaughter->GetLogicalVolume()->TotalVolumeEntities();
  }
  return vols;
}

// ********************************************************************
// GetMass
//
// Returns the mass of the logical volume tree computed from the
// estimated geometrical volume of each solid and material associated
// to the logical volume and its daughters.
// NOTE: the computation may require considerable amount of time,
//       depending from the complexity of the geometry tree.
//       The returned value is cached and can be used for successive
//       calls (default), unless recomputation is forced by providing
//       'true' for the boolean argument in input. Computation should
//       be forced if the geometry setup has changed after the previous
//       call. By setting the 'propagate' boolean flag to 'false' the 
//       method returns the mass of the present logical volume only 
//       (subtracted for the volume occupied by the daughter volumes).
//       The extra argument 'parMaterial' is internally used to
//       consider cases of geometrical parameterisations by material.
// ********************************************************************
//
G4double G4LogicalVolume::GetMass(G4bool forced,
                                  G4bool propagate,
                                  G4Material* parMaterial)
{
  // Return the cached non-zero value, if not forced
  //
  if ( (G4MT_mass) && (!forced) )  { return G4MT_mass; }

  // Global density and computed mass associated to the logical
  // volume without considering its daughters
  //
  G4Material* logMaterial = parMaterial ? parMaterial : GetMaterial();
  if (logMaterial == nullptr)
  {
    std::ostringstream message;
    message << "No material associated to the logical volume: "
            << fName << " !" << G4endl
            << "Sorry, cannot compute the mass ...";
    G4Exception("G4LogicalVolume::GetMass()", "GeomMgt0002",
                FatalException, message);
    return 0.0;
  }
  if ( GetSolid() == nullptr )
  {
    std::ostringstream message;
    message << "No solid is associated to the logical volume: "
            << fName << " !" << G4endl
            << "Sorry, cannot compute the mass ...";
    G4Exception("G4LogicalVolume::GetMass()", "GeomMgt0002",
                FatalException, message);
    return 0.0;
  }
  G4double globalDensity = logMaterial->GetDensity();
  G4double motherMass = GetSolid()->GetCubicVolume() * globalDensity;
  G4double massSum = motherMass;
  
  // For each daughter in the tree, subtract the mass occupied
  // and if required by the propagate flag, add the real daughter's
  // one computed recursively

  for (auto itDau = fDaughters.cbegin(); itDau != fDaughters.cend(); ++itDau)
  {
    G4VPhysicalVolume* physDaughter = (*itDau);
    G4LogicalVolume* logDaughter = physDaughter->GetLogicalVolume();
    G4double subMass = 0.0;
    G4VSolid* daughterSolid = nullptr;
    G4Material* daughterMaterial = nullptr;

    // Compute the mass to subtract and to add for each daughter
    // considering its multiplicity (i.e. replicated or not) and
    // eventually its parameterisation (by solid and/or by material)
    //
    for (auto i=0; i<physDaughter->GetMultiplicity(); ++i)
    {
      G4VPVParameterisation* physParam = physDaughter->GetParameterisation();
      if (physParam)
      {
        daughterSolid = physParam->ComputeSolid(i, physDaughter);
        daughterSolid->ComputeDimensions(physParam, i, physDaughter);
        daughterMaterial = physParam->ComputeMaterial(i, physDaughter);
      }
      else
      {
        daughterSolid = logDaughter->GetSolid();
        daughterMaterial = logDaughter->GetMaterial();
      }
      subMass = daughterSolid->GetCubicVolume() * globalDensity;

      // Subtract the daughter's portion for the mass and, if required,
      // add the real daughter's mass computed recursively
      //
      massSum -= subMass;
      if (propagate)
      {
        massSum += logDaughter->GetMass(true, true, daughterMaterial);
      }
    }
  }
  G4MT_mass = massSum;
  return massSum;
}

// ********************************************************************
// Change the daughters volume type -- checking proposed values
//
//  Undertakes primitive checking, to ensure that only 'legal' changes
//  are made:
//    - any type to 'external'  ( user responsibility )
//    - the type proposed is checked against the deduced type
//       (for potential switch back to 'internal' type.)
//  Returns success (true) or failure (false)
//
G4bool G4LogicalVolume::ChangeDaughtersType(EVolume aType)
{
  G4bool works = false;
  if( aType == kExternal )
  {
    // It is the responsibility of External Navigator to handle types selected
    //
    fDaughtersVolumeType = aType;
    works = true;
  }
  else
  {
    EVolume expectedVType = DeduceDaughtersType();
    works = (expectedVType == aType);
    if ( works )
    {
      fDaughtersVolumeType = aType;
    }
  }
  return works;
}

// ********************************************************************
// SetVisAttributes - copy version
// ********************************************************************
//
void G4LogicalVolume::SetVisAttributes (const G4VisAttributes& VA)
{
  if (G4Threading::IsWorkerThread()) return;
  fVisAttributes = std::make_shared<const G4VisAttributes>(VA);
}

// ********************************************************************
// SetVisAttributes
// ********************************************************************
//
void G4LogicalVolume::SetVisAttributes (const G4VisAttributes* pVA)
{
  if (G4Threading::IsWorkerThread()) return;
  fVisAttributes = std::shared_ptr<const G4VisAttributes>(pVA,[](const G4VisAttributes*){});
}

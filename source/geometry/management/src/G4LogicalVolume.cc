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
// $Id: G4LogicalVolume.cc 90701 2015-06-08 09:43:16Z gcosmo $
//
// 
// class G4LogicalVolume Implementation
//
// History:
// 15.01.13 G.Cosmo, A.Dotti: Modified for thread-safety for MT
// 01.03.05 G.Santin: Added flag for optional propagation of GetMass()
// 17.05.02 G.Cosmo: Added flag for optional optimisation
// 12.02.99 S.Giani: Default initialization of voxelization quality
// 04.08.97 P.M.DeFreitas: Added methods for parameterised simulation 
// 19.08.96 P.Kent: Modified for G4VSensitive Detector
// 11.07.95 P.Kent: Initial version
// --------------------------------------------------------------------

#include "G4LogicalVolume.hh"
#ifdef   NO_INLINE
#define  inline
#include "G4LogicalVolume.icc"
#undef   inline
#endif

#include "G4LogicalVolumeStore.hh"
#include "G4VSolid.hh"
#include "G4Material.hh"
#include "G4VPVParameterisation.hh"
#include "G4VisAttributes.hh"

#include "G4UnitsTable.hh"

// This static member is thread local. For each thread, it points to the
// array of G4LVData instances.
//
template <class G4LVData> G4ThreadLocal
G4LVData* G4GeomSplitter<G4LVData>::offset = 0;      

G4LVData::G4LVData()
: fSolid(0),fSensitiveDetector(0),fFieldManager(0),
  fMaterial(0),fMass(0.),fCutsCouple(0)
{;}

// This new field helps to use the class G4LVManager
//
G4LVManager G4LogicalVolume::subInstanceManager;

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
  AssignFieldManager(fFieldManager); // Should be set - but a per-thread copy is not available yet
  // G4MT_fmanager= fFieldManager;
  //  Must not call SetFieldManager(fFieldManager, false); which propagates FieldMgr

#ifdef CLONE_FIELD_MGR
  // Create a field FieldManager by cloning
  G4FieldManager workerFldMgr= fFieldManager->GetWorkerClone(G4bool* created);
  if( created || (GetFieldManager()!=workerFldMgr) )
  {
    SetFieldManager(fFieldManager, false); // which propagates FieldMgr
  }else{
    // Field manager existed and is equal to current one
    AssignFieldManager(workerFldMgr);
  }
#endif
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
 : fDaughters(0,(G4VPhysicalVolume*)0), 
   fVoxel(0), fOptimise(optimise), fRootRegion(false), fLock(false),
   fSmartless(2.), fVisAttributes(0), fRegion(0), fBiasWeight(1.)
{
  // Initialize 'Shadow'/master pointers - for use in copying to workers
  fSolid = pSolid;
  fSensitiveDetector = pSDetector;
  fFieldManager = pFieldMgr;

  instanceID = subInstanceManager.CreateSubInstance();
  AssignFieldManager(pFieldMgr); // G4MT_fmanager = pFieldMgr;
  
  // fMasterFieldMgr= pFieldMgr;
  G4MT_mass = 0.;
  G4MT_ccouple = 0;

  SetSolid(pSolid);
  SetMaterial(pMaterial);
  SetName(name);
  SetSensitiveDetector(pSDetector);
  SetUserLimits(pULimits);    
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
 : fDaughters(0,(G4VPhysicalVolume*)0),
   fName(""), fUserLimits(0),
   fVoxel(0), fOptimise(true), fRootRegion(false), fLock(false),
   fSmartless(2.), fVisAttributes(0), fRegion(0), fBiasWeight(1.),
   fSolid(0), fSensitiveDetector(0), fFieldManager(0)
{
  instanceID = subInstanceManager.CreateSubInstance();
  
  SetSensitiveDetector(0);    // G4MT_sdetector = 0;
  SetFieldManager(0, false);  // G4MT_fmanager = 0;

  G4MT_mass = 0.;
  G4MT_ccouple = 0;
  
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
  G4LogicalVolumeStore::DeRegister(this);
}

// ********************************************************************
// SetFieldManager
// ********************************************************************
//
void
G4LogicalVolume::SetFieldManager(G4FieldManager* pNewFieldMgr,
                                 G4bool          forceAllDaughters) 
{
  // G4MT_fmanager = pNewFieldMgr;
  AssignFieldManager(pNewFieldMgr);

  G4int NoDaughters = GetNoDaughters();
  while ( (NoDaughters--)>0 )
  {
    G4LogicalVolume* DaughterLogVol; 
    DaughterLogVol = GetDaughter(NoDaughters)->GetLogicalVolume();
    if ( forceAllDaughters || (DaughterLogVol->GetFieldManager() == 0) )
    {
      DaughterLogVol->SetFieldManager(pNewFieldMgr, forceAllDaughters);
    }
  }
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
    for (G4PhysicalVolumeList::const_iterator itDau = fDaughters.begin();
         itDau != fDaughters.end(); itDau++)
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
  for (G4PhysicalVolumeList::const_iterator itDau = fDaughters.begin();
       itDau != fDaughters.end(); itDau++)
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
  if ( (G4MT_mass) && (!forced) ) return G4MT_mass;

  // Global density and computed mass associated to the logical
  // volume without considering its daughters
  //
  G4Material* logMaterial = parMaterial ? parMaterial : GetMaterial(); // G4MT_material;
  if (!logMaterial)
  {
    std::ostringstream message;
    message << "No material associated to the logical volume: " << fName << " !"
            << G4endl
            << "Sorry, cannot compute the mass ...";
    G4Exception("G4LogicalVolume::GetMass()", "GeomMgt0002",
                FatalException, message);
    return 0;
  }
  if (! GetSolid() ) // !G4MT_solid)
  {
    std::ostringstream message;
    message << "No solid is associated to the logical volume: " << fName << " !"
            << G4endl
            << "Sorry, cannot compute the mass ...";
    G4Exception("G4LogicalVolume::GetMass()", "GeomMgt0002",
                FatalException, message);
    return 0;
  }
  G4double globalDensity = logMaterial->GetDensity();
  G4double motherMass= GetSolid()->GetCubicVolume() * globalDensity;

  // G4MT_mass =
  //  SetMass( motherMmass );
  G4double massSum= motherMass;
  
  // For each daughter in the tree, subtract the mass occupied
  // and if required by the propagate flag, add the real daughter's
  // one computed recursively

  for (G4PhysicalVolumeList::const_iterator itDau = fDaughters.begin();
       itDau != fDaughters.end(); itDau++)
  {
    G4VPhysicalVolume* physDaughter = (*itDau);
    G4LogicalVolume* logDaughter = physDaughter->GetLogicalVolume();
    G4double subMass=0.;
    G4VSolid* daughterSolid = 0;
    G4Material* daughterMaterial = 0;

    // Compute the mass to subtract and to add for each daughter
    // considering its multiplicity (i.e. replicated or not) and
    // eventually its parameterisation (by solid and/or by material)
    //
    for (G4int i=0; i<physDaughter->GetMultiplicity(); i++)
    {
      G4VPVParameterisation*
        physParam = physDaughter->GetParameterisation();
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
  G4MT_mass= massSum;
  return massSum;
}

void G4LogicalVolume::SetVisAttributes (const G4VisAttributes& VA)
{
  fVisAttributes = new G4VisAttributes(VA);
}

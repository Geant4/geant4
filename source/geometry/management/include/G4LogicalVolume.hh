// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LogicalVolume.hh,v 1.3 1999-11-11 15:35:49 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4LogicalVolume
//
// Represents a leaf node or unpositioned subtree in the geometry hierarchy.
// Logical volumes are named, and may have daughters ascribed to them.
// They are responsible for retrieval of the physical and tracking attributes
// of the physical volume that it represents: Solid, material, magnetic field,
// and optionally: user limits, sensitive detectors.
//
// Get and Set functionality is provided for all atributes, but note that
// most set functions should not be used  when the geometry is `closed'.
// As a  further development, `Guard' checks can be added to ensure
// only legal operations at tracking time.
//
// On construction, solid, material and name must be specified
//
//
// Daughters are ascribed and managed by means of a simple
// GetNoDaughters,Get&SetDaughter(n),AddDaughter interface
//
// Smart voxels as used for tracking optimisation are also an attribute.
//
// Logical volumes self register to the logical volume Store on construction,
// and deregister on destruction.
//
// NOTE: This class is currently *NOT* subclassed. If subclassed make
// destructor virtual.
//
// Member functions:
//
// G4LogicalVolume(const G4VSolid *pSolid, const G4Material *pMaterial,
//                 const G4String& name,
//                 const G4MagneticField *pField=0,
//                 const G4VSensitiveDetector *pSDetector=0,
//                 const G4UserLimits *pULimits=0)
//
//   Constructor. The solid and material pointer must be non null. The
//   parameters for field, detector and user limits are optional.
//   The volume also enters itself into the logical volume Store.
//
// ~G4LogicalVolume()
//   Destructor. Removes the logical volume from the logical volume Store.
//
// G4String GetName() const
//   Returns name of logical volume
// void SetName(const G4String& pName)
//   Sets name of logical volume
//
// G4int GetNoDaughters() const
//   Returns the number of daughters (0 to n)
// G4VPhysicalVolume* GetDaughter(const G4int i) const
//   Return the ith daughter. Note numbering starts from 0, and no bounds
//   checkingis performed.
// void SetDaughter(const G4int i,G4VPhysicalVolume* p) 
//   Set the ith daughter to be p, where 0<=i<GetNoDaughters(). Intended
//   for UI use only
// void AddDaughter(G4VPhysicalVolume* p)
//   Add the volume p as a daughter of the current logical volume.
// G4bool IsDaughter(const G4VPhsyicalVolume* p) const
//   Returns true is the volume p is a daughter of the current logical volume
// void RemoveDaughter(const G4VPhysicalVolume* p )
//   Remove the volume p from the List of daughter of the current logical
//   volume.
//
// G4VSolid* GetSolid() const
//   Gets current solid.
// void SetSolid(G4VSolid *pSolid)
//   Sets solid.
//
// G4Material* GetMaterial() const
//   Gets current Material.
// void SetMaterial(G4Material *pMaterial)
//   Sets Material.
//
// G4FieldManager* GetFieldManager() const
//   Gets current FieldManager.
// void SetFieldManager(G4FieldManager *pField, G4bool forceToAllDaughters)
//   Sets FieldManager and propagates it
//    i) only to daughters with G4FieldManager = 0 if forceToAllDaughters=false
//   ii) to all daughters if forceToAllDaughters=true
//
// G4VSensitiveDetector* GetSensitiveDetector() const
//   Gets current SensitiveDetector.
// void SetSensitiveDetector(G4VSensitiveDetector *pSDetector)
//   Sets SensitiveDetector (can be NULL)
//
// G4UserLimits* GetUserLimits() const
//   Gets current UserLimits.
// void SetUserLimits(G4UserLimits *pULimits)
//   Sets UserLimits.
//
// G4VoxelHeader* GetVoxelHeader() const
//   Gets current VoxelHeader.
// void SetVoxelHeader(G4VoxelHeader *pVoxel)
//   Sets VoxelHeader.
//
// G4double GetSmartless()
//   Gets user defined optimisation quality
// void SetSmartless(G4double)
//   Sets user defined optimisation quality
//
// void BecomeEnvelopeForFastSimulation(G4FastSimulationManager* );
//   Makes this an Envelope for given FastSimulationManager. 
//   Ensures that all its daughter volumes get it too - unless they 
//   have one already.
// G4FastSimulationManager* GetFastSimulationManager () const;
//   Gets current FastSimulationManager pointer.
// void ClearEnvelopeForFastSimulation(G4LogicalVolume* motherLogVol);
//   Erase volume's Envelope status and propagate the FastSimulationManager 
//   of its mother volume to itself and its daughters.
//
// void  SetFastSimulationManager (G4FastSimulationManager* pPA, 
//                           G4bool IsEnvelope);
//   Sets the fast simulation manager. Private method called by the
//   public SetIsEnvelope method with IsEnvelope = TRUE. It is 
//   then called recursivaly to the daughters to propagate the 
//   FastSimulationManager pointer with IsEnvelope = FALSE.
//
// void SetBiasWeight (G4double w);
//   Sets the bias weight
// G4double GetBiasWeight() const;
//   Gets the bias weight
//
// Operators:
//
// G4bool operator == (G4LogicalVolume,G4LogicalVolume)
//   Equality defined by address only- return true if objects are at
//   same address, else false
//
//
// Member data:
//
// G4RWTPtrOrderedVector<G4VPhysicalVolume> fDaughters
//   Vector of daughters. Given initial size of 0.
// G4FieldManager *fFieldManager 
//   Pointer (possibly NULL) to (magnetic or other) field manager object
// G4Material *fMaterial
//   Pointer to material at this node
// G4String fName
//   Name of logical volume
// G4SensitiveDetector *fSensitiveDetector
//   Pointer (possibly NULL) to `Hit' object
// G4VSolid *fSolid
//   Pointer to solid
// G4UserLimits *fUserLimits
//   Pointer (possibly NULL) to user Step limit object for this node
// G4VoxelHeader *fVoxel
//   Pointer (possibly NULL) to optimisation info objects
// G4double smartless
//   Quality for optimisation, average number of voxels to be spent per content 
// G4FastSimulationManager *fFastSimulationManager
//   Pointer (possibly NULL) to G4FastSimulationManager object
// G4bool fIsEnvelope
//   Flags if the Logical Volume is an envelope for a FastSimulationManager.
// G4double fBiasWeight
//   weight used in the event biasing technique
//
// History:
// 12.02.99 S.Giani: Added user defined optimisation quality
// 09.11.98 J. Apostolakis:  Changed G4MagneticField to G4FieldManager
// 09.11.98 M. Verderi & JA. : added  BiasWeight member and Get/Set methods
// 10.20.97 P. MoraDeFreitas : added pointer to a FastSimulation
//          (J.Apostolakis)      & flag to indicate if it is an Envelope for it
// 19.11.96 J.Allison Replaced G4Visible with explicit const G4VisAttributes*.
// 19.08.96 P.Kent    Split -> hh/icc/cc files; G4VSensitiveDetector change
// 11.07.95 P.Kent    Initial version.

#ifndef G4LOGICALVOLUME_HH
#define G4LOGICALVOLUME_HH

#include "globals.hh"
#include "G4VPhysicalVolume.hh"	// Need operator == for vector fdaughters
#include "g4rw/tpordvec.h"
#include <assert.h>

// Forward declarations
class G4FieldManager;
class G4Material;
class G4VSensitiveDetector;
class G4VSolid;
class G4UserLimits;
class G4SmartVoxelHeader;
class G4VisAttributes;
class G4FastSimulationManager;

class G4LogicalVolume
{
public:
    G4LogicalVolume(G4VSolid *pSolid, G4Material *pMaterial,
		    const G4String& name,
		    G4FieldManager *pFieldMgr=0,
		    G4VSensitiveDetector *pSDetector=0,
		    G4UserLimits *pULimits=0);

    ~G4LogicalVolume();

    G4String GetName() const;
    void SetName(const G4String& pName);

    G4int GetNoDaughters() const;
    G4VPhysicalVolume* GetDaughter(const G4int i) const;
    void AddDaughter(G4VPhysicalVolume* p);
    G4bool IsDaughter(const G4VPhysicalVolume* p) const;
    void RemoveDaughter(const G4VPhysicalVolume* p);

    G4VSolid* GetSolid() const;
    void SetSolid(G4VSolid *pSolid);
    
    G4Material* GetMaterial() const;
    void SetMaterial(G4Material *pMaterial);

    G4FieldManager* GetFieldManager() const;
    void SetFieldManager(G4FieldManager *pFieldMgr, G4bool forceToAllDaughters); 

    G4VSensitiveDetector* GetSensitiveDetector() const;
    void SetSensitiveDetector(G4VSensitiveDetector *pSDetector);

    G4UserLimits* GetUserLimits() const;
    void SetUserLimits(G4UserLimits *pULimits);

    G4SmartVoxelHeader* GetVoxelHeader() const;
    void SetVoxelHeader(G4SmartVoxelHeader *pVoxel);
    
    G4double GetSmartless();
    void SetSmartless(G4double s);
    
    G4bool operator == ( const G4LogicalVolume& lv) const;

    const G4VisAttributes* GetVisAttributes () const;
    void  SetVisAttributes (const G4VisAttributes* pVA);
    void  SetVisAttributes (const G4VisAttributes& VA);

    void BecomeEnvelopeForFastSimulation(G4FastSimulationManager* );
    void  ClearEnvelopeForFastSimulation(G4LogicalVolume* motherLV= 0);
    G4FastSimulationManager* GetFastSimulationManager () const;

    void SetBiasWeight (G4double w);
    G4double GetBiasWeight() const;

private:
    void  SetFastSimulationManager (G4FastSimulationManager* pPA, 
			      G4bool IsEnvelope);
    G4LogicalVolume* FindMotherLogicalVolumeForEnvelope(); 
//   
//  Data members:   

private:
    G4RWTPtrOrderedVector<G4VPhysicalVolume> fDaughters;    
    G4FieldManager *fFieldManager;
    G4Material *fMaterial;
    G4String fName;
    G4VSensitiveDetector *fSensitiveDetector;
    G4VSolid *fSolid;
    G4UserLimits *fUserLimits;
    G4SmartVoxelHeader *fVoxel;
    G4double fSmartless;
    const G4VisAttributes* fVisAttributes;
    G4FastSimulationManager *fFastSimulationManager;
    G4bool fIsEnvelope;
    G4double fBiasWeight;
};

#include "G4LogicalVolume.icc"

#endif

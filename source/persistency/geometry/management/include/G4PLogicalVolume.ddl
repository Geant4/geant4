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
// $Id: G4PLogicalVolume.ddl,v 1.11 2001/07/11 10:02:18 gunter Exp $
// GEANT4 tag $Name: geant4-04-00 $

// Class Description:
//   Persistent version of the class G4LogicalVolume.

#ifndef G4PLOGICALVOLUME_DDL
#define G4PLOGICALVOLUME_DDL 1

#include "G4Pglobals.hh"
#include "G4PersistentTypes.hh" 
#include "G4PersistentSchema.hh"
#include "G4VMaterialMap.hh"
#include "HepODBMS/odbms/HepODBMS.h"

#include "G4PVSolid.hh"

// Forward declarations
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSolid;

#include "G4PVPhysicalVolume.hh"

typedef d_Varray< d_Ref<G4PVPhysicalVolume> > G4PVPhysVolRefVArray;


class G4PLogicalVolume: public HepPersObj
{
  public: // With description
    G4PLogicalVolume( const G4LogicalVolume* aLogicalVolume,
                      HepRef(G4PVSolid) persSolid);
    ~G4PLogicalVolume();
      // Constructor and Destructor

    G4LogicalVolume* MakeTransientObject(G4VSolid* aSolid,
                                         G4Material* aMaterial);
      // Creates a transient G4LogicalVolume object.

    HepRef(G4PVSolid) GetSolid();
    void SetSolid( HepRef(G4PVSolid) pSolid );
      // Get or set the pointer of persistent solid object.

  public: // With description
    G4int GetNoDaughters() const;
    HepRef(G4PVPhysicalVolume) GetDaughter(const G4int i) const;
    void AddDaughter(HepRef(G4PVPhysicalVolume) p);
    G4bool IsDaughter(const HepRef(G4PVPhysicalVolume) p) const;
    void RemoveDaughter(const HepRef(G4PVPhysicalVolume) p);
      // Methods to handle association of persistent physical volume.

  private:
    G4PVPhysVolRefVArray fDaughters;    
//    G4MagneticField* fMagneticField;
//    G4Material* fMaterial;
    G4PString fMaterialName;
    G4PString fName;
//    G4VSensitiveDetector *fSensitiveDetector;
    d_Ref<G4PVSolid>    fSolid;
//    G4UserLimits* fUserLimits;
//    G4SmartVoxelHeader* fVoxel;
//    const G4VisAttributes* fVisAttributes;
//    G4FastSimulationManager* fFastSimulationManager;
//    G4bool fIsEnvelope;

    G4Pint NoOfDaughters;

  public: // With description
      inline void SetMaterialName(G4String aName)
      { fMaterialName = aName; }
      inline G4String GetMaterialName()
      { return (G4String) fMaterialName; }
        // Set and Get methods for the associated material name.

};

#endif

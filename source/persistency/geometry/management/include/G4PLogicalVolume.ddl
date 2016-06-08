// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PLogicalVolume.ddl,v 1.8 2000/06/09 12:56:01 morita Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
// 
// class G4PLogicalVolume
// P-verson of the class G4LogicalVolume
//                                              Takashi.Sasaki@kek.jp
#ifndef G4PLOGICALVOLUME_DDL
#define G4PLOGICALVOLUME_DDL 1

#include "globals.hh"
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
public:
    G4PLogicalVolume( const G4LogicalVolume* aLogicalVolume,
                      HepRef(G4PVSolid) persSolid);
    ~G4PLogicalVolume();

    G4LogicalVolume* MakeTransientObject(G4VSolid* aSolid,
                                         G4Material* aMaterial);

    HepRef(G4PVSolid) GetSolid();
    void SetSolid( HepRef(G4PVSolid) pSolid );

public:
    G4int GetNoDaughters() const;
    HepRef(G4PVPhysicalVolume) GetDaughter(const G4int i) const;
    void AddDaughter(HepRef(G4PVPhysicalVolume) p);
    G4bool IsDaughter(const HepRef(G4PVPhysicalVolume) p) const;
    void RemoveDaughter(const HepRef(G4PVPhysicalVolume) p);

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

  public:
      inline void SetMaterialName(G4String aName)
      { fMaterialName = aName; }
      inline G4String GetMaterialName()
      { return (G4String) fMaterialName; }
};

#endif

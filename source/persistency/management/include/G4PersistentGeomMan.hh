// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PersistentGeomMan.hh,v 1.1 1999/01/07 16:10:56 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
// class G4PersistentGeomMan 
//
// A Utility class for storing and retrieving the geometry objects.
//
// This class is not persistent-capable. 
//
// Member functions:
// =================
//
//
// Member data:
// ============
// static G4PersistentGeomMan* fPersistencyManager
//   Ptr to the unique instance of class
//
// History:
// 98.10.30 Y.Morita  Splited from G4PersistencyManager

#ifndef G4PERSISTENTGEOMMAN_H
#define G4PERSISTENTGEOMMAN_H 1

#include "globals.hh"
#include "geomdefs.hh"

#include "HepODBMS/clustering/HepDbApplication.h"

#define G4_PHYS_VOLUME_DEPTH_MAX 1000000

class G4VSolid;
class G4PEvent;
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4PVPhysicalVolume;
class G4PLogicalVolume;
class G4PVSolid;

#include "G4PGeometryObjectMap.hh"


class G4PersistentGeomMan 
{
  public:
      G4PersistentGeomMan();
      ~G4PersistentGeomMan();

  private:
      HepDatabaseRef f_GeomDB;
      HepContainerRef f_GeomContainer;
      HepRef(G4PGeometryObjectMap) f_GeomMap;

  public:
      G4bool Store( HepDbApplication* dbApp,
                    const G4VPhysicalVolume* aWorld );
      G4bool Retrieve( HepDbApplication* dbApp,
                       G4VPhysicalVolume*& aWorld );

  private:
      // Characterise `type' of volume - normal/replicated/parameterised
      EVolume VolumeType(const G4VPhysicalVolume *pVol) const;

      HepRef(G4PVPhysicalVolume) MakePersistentObject (
                                        G4VPhysicalVolume* aPhysVol );
      HepRef(G4PLogicalVolume) MakePersistentObject (
                                        G4LogicalVolume* aLogVol );
      HepRef(G4PVSolid) MakePersistentObject ( G4VSolid* aSolid );

      G4VPhysicalVolume* MakeTransientObject (
                             HepRef(G4PVPhysicalVolume) aPhysVol,
                             HepRef(G4PVPhysicalVolume) theMotherVol );
      G4LogicalVolume* MakeTransientObject (
                             HepRef(G4PLogicalVolume) aLogVol,
                             HepRef(G4PVPhysicalVolume) theMotherVol );
      G4VSolid* MakeTransientObject ( HepRef(G4PVSolid) aSolid );

  private:
      G4int f_nRecursive;
      G4int f_tRecursive;
      G4PString f_GeomMapScopeName;
};

inline EVolume G4PersistentGeomMan::VolumeType(const G4VPhysicalVolume *pVol) const
{
  EVolume type;
  EAxis axis;
  G4int nReplicas;
  G4double width,offset;
  G4bool consuming;
  if (pVol->IsReplicated())
    {
      pVol->GetReplicationData(axis,nReplicas,width,offset,consuming);
      type=(consuming) ? kReplica : kParameterised;
    }
  else
    {
      type=kNormal;
    }
  return type;
}

#endif


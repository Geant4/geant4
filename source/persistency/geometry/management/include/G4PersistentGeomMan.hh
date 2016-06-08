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
// $Id: G4PersistentGeomMan.hh,v 1.11 2001/07/11 10:02:19 gunter Exp $
// GEANT4 tag $Name: geant4-04-01 $
//

// Class Description:
//   A Utility class to be used by G4PersistencyManager.
// Average users do not need to use this class.
//
// This class is not persistent-capable.
//

// History:
// 98.10.30 Y.Morita  Splited from G4PersistencyManager

#ifndef G4PERSISTENTGEOMMAN_H
#define G4PERSISTENTGEOMMAN_H 1

#include "globals.hh"
#include "geomdefs.hh"

#include "HepODBMS/odbms/HepODBMS.h"

#include "G4VPersistentSubMan.hh"
#include "G4VPersistentSubDbMan.hh"
#include "G4PGeometryObjectMap.hh"
#include "G4VMaterialMap.hh"

#define G4_PHYS_VOLUME_DEPTH_MAX 1000000

// forward declarations
class HepDbApplication;
class G4VSolid;
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4PVPhysicalVolume;
class G4PLogicalVolume;
class G4PVSolid;


class G4PersistentGeomMan 
 : public G4VPersistentSubMan, public G4VPersistentSubDbMan 
{
  friend class G4PersistencyManager;
  friend class G4TransactionManager;

  private:
      // to be used by G4PersistencyManager only 
      G4PersistentGeomMan();
      ~G4PersistentGeomMan();

  private:
      HepRef(G4PGeometryObjectMap) f_GeomMap;
      G4VMaterialMap* f_MaterialMap;

  private:
      // interface with G4PersistencyManager
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
      G4String f_GeomMapScopeName;
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


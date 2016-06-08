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
// $Id: G4PGeometryObjectMap.ddl,v 1.9 2001/07/11 10:02:18 gunter Exp $
// GEANT4 tag $Name: geant4-04-00 $
//

// class description:
//
//	This is a class which is responsible for keeping track of the 
//	persistent geometry object in storing and retrieving geometry.
//	This class inherits HepPersObj and is persistent-capable. 
//

// Note:
// =====
//      This class contains pointers of some transient classes.
//      As a result, ooddlx compiler will issue warnings as follows.
//      The value of the pointers are used only in transient case
//      in this class, and you can safely ignore these warnings.
// ----
// "include/G4PGeometryObjectMap.ddl", line 136: warning: persistent-capable
//           class member type is pointer
//         G4VPhysVolRefArray*  transPhysVolPtrs;
//                              ^
// ----
//
// History:
// 98.06.20 Y.Morita  Initial version

#ifndef G4PGeometryObjectMap_h
#define G4PGeometryObjectMap_h 1

#include "G4Pglobals.hh"

#include "g4std/vector"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4PersistentTypes.hh"
#include "G4PersistentSchema.hh"
#include "HepODBMS/odbms/HepODBMS.h"

class G4PVPhysicalVolume;
#pragma ooclassref G4PVPhysicalVolume "G4PVPhysicalVolume_ref.hh"
class G4PLogicalVolume;
#pragma ooclassref G4PLogicalVolume "G4PLogicalVolume_ref.hh"
class G4PVSolid;
#pragma ooclassref G4PVSolid "G4PVSolid_ref.hh"

typedef G4std::vector<G4LogicalVolume*>       G4LogVolRefArray;
typedef G4std::vector<G4VPhysicalVolume*>     G4VPhysVolRefArray;
typedef G4std::vector<G4VSolid*>              G4VSolidRefArray;

typedef G4std::vector<G4LogicalVolume*>::iterator   G4LogVolRefArrayItr;
typedef G4std::vector<G4VPhysicalVolume*>::iterator G4VPhysVolRefArrayItr;
typedef G4std::vector<G4VSolid*>::iterator          G4VSolidRefArrayItr;

typedef d_Varray< d_Ref<G4PVPhysicalVolume> > G4PVPhysVolRefVArray;
typedef d_Varray< d_Ref<G4PLogicalVolume> >   G4PLogVolRefVArray;
typedef d_Varray< d_Ref<G4PVSolid> >          G4PVSolidRefVArray;

class G4PGeometryObjectMap 
 : public HepPersObj
{
  public:
      G4PGeometryObjectMap();
  public: // with description
      G4PGeometryObjectMap( const G4String theGeometryName );
      ~G4PGeometryObjectMap();
      //  The constructor and the destructor.

  public: // with description
      HepRef(G4PVPhysicalVolume) LookUp(         G4VPhysicalVolume* aPhysVol );
      G4VPhysicalVolume*         LookUp( HepRef(G4PVPhysicalVolume) aPhysVol );
      HepRef(G4PLogicalVolume)   LookUp(         G4LogicalVolume*   aLogVol);
      G4LogicalVolume*           LookUp( HepRef(G4PLogicalVolume)   aLogVol);
      HepRef(G4PVSolid)          LookUp(         G4VSolid*          aSolid);
      G4VSolid*                  LookUp( HepRef(G4PVSolid)          persSolid);
      // Method LookUp() will check if the pair of transient/persistent
      // geometry objects is already registered in the geometry map.
      // It returns the pointer of the pairing object if the pair exists.

  public: // with description
      void Add( G4VPhysicalVolume*         aPhysVol,
                HepRef(G4PVPhysicalVolume) persPhysVol );
      void Add( HepRef(G4PVPhysicalVolume) persPhysVol,
                G4VPhysicalVolume*         aPhysVol );
      void Add( G4LogicalVolume*           aLogVol,
                HepRef(G4PLogicalVolume)   persLogVol );
      void Add( HepRef(G4PLogicalVolume)   persLogVol,
                G4LogicalVolume*           aLogVol );
      void Add( G4VSolid*                  aSolid,
                HepRef(G4PVSolid)          persSolid );
      void Add( HepRef(G4PVSolid)          persSolid,
                G4VSolid*                  aSolid );
      // Method Add() will register the transient and persistent geometry
      // objects into the geometry map.

  public: // with description
      inline void SetWorldVolume( const HepRef(G4PVPhysicalVolume) aWorld )
            { thePersistentWorldVolume = aWorld; }
      inline HepRef(G4PVPhysicalVolume) GetWorldVolume()
            { return thePersistentWorldVolume; }
      // Set and Get methods for setting and getting a smart pointer
      // for the persistent world volume of the geometry.

  public: // with description
      inline G4int GetNoPhysVol() const { return noPhysVol; }
      inline G4int GetNoLogVol()  const { return noLogVol; }
      inline G4int GetNoSolids()  const { return noSolids; }
      // returns the number of the registered transient/persistent pairs
      // in the geometry map.

  public: // with description
      G4VSolid*         GetSolid(G4int n);
      HepRef(G4PVSolid) GetPSolid(G4int n);
      // returns the pointer of n-th transient/persistent solid.

  public:
      void InitTransientMap();

  private:
      // World Volume of this Geometry object tree
      d_Ref<G4PVPhysicalVolume> thePersistentWorldVolume;

      // Physics Volume Pointer Array
      G4Pint noPhysVol;
      G4VPhysVolRefArray*  transPhysVolPtrs;
      G4PVPhysVolRefVArray persPhysVolPtrs;

      // Logical Volume Pointer Array
      G4Pint noLogVol;
      G4LogVolRefArray*  transLogVolPtrs;
      G4PLogVolRefVArray persLogVolPtrs;

      // Solid Pointer Array
      G4Pint noSolids;
      G4VSolidRefArray*  transSolidPtrs;
      G4PVSolidRefVArray persSolidPtrs;

};

#endif

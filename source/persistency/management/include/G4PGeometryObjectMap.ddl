// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PGeometryObjectMap.ddl,v 2.5 1998/07/14 14:12:35 morita Exp $
// GEANT4 tag $Name: geant4-00 $
//
// class G4PGeometryObjectMap 
//
// A Class responsible for keeping track of the geometry object map.
//
// This class is persistent-capable. 
//
// Member functions:
// =================
//      G4PGeometryObjectMap( G4PString theGeometryName );
//      ~G4PGeometryObjectMap();
//
//  protected:
//      G4PVPhysicalVolume* LookUp( G4VPhysicalVolume* aPhysVol );
//      void Add( const G4VPhysicalVolume*  aPhysVol,
//                const G4PVPhysicalVolume* persPhysVol );
//      G4PLogicalVolume* LookUp(G4LogicalVolume* aLogVol);
//      void Add( const G4LogicalVolume*  aLogVol,
//                const G4PLogicalVolume* persLogVol );
//      G4PVSolid* LookUp(G4VSolid* aSolid);
//      void Add( const G4VSolid*  aSolid,
//                const G4PVSolid* persSolid );
//      G4int GetNoSolids();
//      G4VSolid* GetSolid(G4int n);
//      G4PVSolid* GetPSolid(G4int n);
//
// Note:
// =====
//      This class contains RWTPtrVector of the transient classes.
//      As a result, ooddlx compiler will issue warnings as follows, but it's okay.
//      The value of RWTPtrVector is used only in transient case
//      in this class, and you can safely ignore these warnings.
// ----
// "include/G4PGeometryObjectMap.ddl", line xxx: warning: persistent-capable
//           class member type is a problem
//             The type of base class RWPtrVector is a problem
//             The type of member RWPtrVector::array_ is pointer
//             The type of base class RWTPtrVector<G4VPhysicalVolume > is a
//                       problem
//             The type of base class RWPtrVector is a problem
//         G4PhysVolRefVArray  transPhysVolPtrs;
//                             ^
// ----
//
//
// History:
// 98.06.20 Y.Morita  Initial version

#ifndef G4PGeometryObjectMap_h
#define G4PGeometryObjectMap_h 1

#include <rw/tpordvec.h>

#include "globals.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4PersistentTypes.hh"
#include "HepODBMS/odbms/HepODBMS.h"

#include "G4PVPhysicalVolume.hh"
#include "G4PVPhysVolRefVArray.hh"
#include "G4PLogicalVolume.hh"
#include "G4PLogVolRefVArray.hh"
#include "G4PVSolid.hh"
#include "G4PVSolidRefVArray.hh"

typedef RWTPtrVector<G4VPhysicalVolume> G4VPhysVolRefVArray;
typedef RWTPtrVector<G4LogicalVolume> G4LogVolRefVArray;
typedef RWTPtrVector<G4VSolid> G4VSolidRefVArray;

class G4PGeometryObjectMap 
 : public HepPersObj
{
  public:
      G4PGeometryObjectMap();
      G4PGeometryObjectMap( const G4String theGeometryName );
      ~G4PGeometryObjectMap();

//    template cannot be used in class scope!!
//      template<class T_IN, class T_OUT> T_OUT LookUp( const T_IN aGeomObj );

      HepRef(G4PVPhysicalVolume) LookUp(         G4VPhysicalVolume* aPhysVol );
      G4VPhysicalVolume*         LookUp( HepRef(G4PVPhysicalVolume) aPhysVol );
      HepRef(G4PLogicalVolume)   LookUp(         G4LogicalVolume*   aLogVol);
      G4LogicalVolume*           LookUp( HepRef(G4PLogicalVolume)   aLogVol);
      HepRef(G4PVSolid)          LookUp(         G4VSolid*          aSolid);
      G4VSolid*                  LookUp( HepRef(G4PVSolid)          persSolid);

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

      inline void SetWorldVolume( const HepRef(G4PVPhysicalVolume) aWorld )
            { thePersistentWorldVolume = aWorld; }
      inline HepRef(G4PVPhysicalVolume) GetWorldVolume()
            { return thePersistentWorldVolume; }

      inline G4int GetNoPhysVol() const { return noPhysVol; }
      inline G4int GetNoLogVol()  const { return noLogVol; }
      inline G4int GetNoSolids()  const { return noSolids; }

      G4VSolid*         GetSolid(G4int n);
      HepRef(G4PVSolid) GetPSolid(G4int n);

      void InitTransientMap();

  private:
      // World Volume of this Geometry object tree
      ooRef(G4PVPhysicalVolume) thePersistentWorldVolume;

      // Physics Volume Pointer Array
      G4Pint noPhysVol;
      G4VPhysVolRefVArray  transPhysVolPtrs;
      G4PVPhysVolRefVArray persPhysVolPtrs;

      // Logical Volume Pointer Array
      G4Pint noLogVol;
      G4LogVolRefVArray  transLogVolPtrs;
      G4PLogVolRefVArray persLogVolPtrs;

      // Solid Pointer Array
      G4Pint noSolids;
      G4VSolidRefVArray  transSolidPtrs;
      G4PVSolidRefVArray persSolidPtrs;

};

#endif

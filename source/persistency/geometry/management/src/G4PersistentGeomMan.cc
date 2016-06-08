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
// $Id: G4PersistentGeomMan.cc,v 1.15.2.1 2001/06/28 19:11:29 gunter Exp $
// GEANT4 tag $Name:  $
//
// class G4PersistentGeomMan 
//
// Implementation for concrete G4PersistentGeomMan.
//
// History:
// 98.10.30 Y.Morita  Splited from G4PersistencyManager

#include <assert.h>
#include "G4PersistentGeomMan.hh"

// include files for forward declarations
#include "HepODBMS/clustering/HepDbApplication.h"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPhysicalVolume.hh"
#include "G4PLogicalVolume.hh"
#include "G4PVSolid.hh"
#include "G4PGeometryObjectMap.hh"

// other used classes
#include "G4PBox.hh"
#include "G4PCons.hh"
#include "G4PPara.hh"
#include "G4PSphere.hh"
#include "G4PTorus.hh"
#include "G4PTrap.hh"
#include "G4PTrd.hh"
#include "G4PTubs.hh"

//#include "G4PEllipticalTube.hh"
#include "G4PHype.hh"
//#include "G4PPolycone.hh"
//#include "G4PPolyhedra.hh"

#include "G4PDisplacedSolid.hh"
#include "G4PIntersectionSolid.hh"
#include "G4PSubtractionSolid.hh"
#include "G4PUnionSolid.hh"

#include "G4BooleanSolid.hh"
#include "G4DisplacedSolid.hh"

#include "G4PPVPlacement.hh"
#include "G4PPVReplica.hh"
#include "G4PPVParameterised.hh"

#include "G4ios.hh"

G4PersistentGeomMan::G4PersistentGeomMan()
 : f_GeomMap(0), f_GeomMapScopeName("Geant4 Geometry Object Map")
{;}

G4PersistentGeomMan::~G4PersistentGeomMan()
{;}

//----------------------------------------------------------------------------

G4bool G4PersistentGeomMan::Store( HepDbApplication* dbApp,
                                   const G4VPhysicalVolume* theWorld)
{
  void* aWorld = (void*) theWorld;

  const G4String theGeometryName = "Geant4 geometry";

  // clear the recursive call counter
  f_nRecursive = 0;
  f_tRecursive = 0;

  // ----------------------- Objectivity Feature ----------------------- //
  // Look for the geometry object map in the database and create it
  // if it does not exist
  if ( ! f_GeomMap.lookupObj( f_DB, f_GeomMapScopeName, oocUpdate ) )
  {

    // Prepare lookup table for the geometry object tree
    f_GeomMap =
       new(f_container) G4PGeometryObjectMap( theGeometryName );

    // Add ScopeName to f_GeomMap for future reference (*Objy feature)
    if ( ! f_GeomMap.nameObj( f_DB, f_GeomMapScopeName ) )
    {
      G4cerr << "G4PersistentGeomMan::Store"
             << " -- error in naming the geometry object map: "
             << f_GeomMapScopeName << G4endl;
      return false;
    }

  }
  else
  {
    G4cerr << "G4PersistentGeomMan: Geometry alreday stored in the"
           << "geometry database." << G4endl
           << "  Geometry not stored." << G4endl;
    return false;
  }
  // ----------------------- Objectivity Feature ----------------------- //

  // make the persistent world volume object from "theWorld"
  HepRef(G4PVPhysicalVolume) persWorld =
                  MakePersistentObject( (G4VPhysicalVolume*)aWorld );

  if ( persWorld != 0 )
  {
    // set the persistent World Volume to f_GeomMap
    f_GeomMap->SetWorldVolume( persWorld );
    return true;
  }
  else
  {
    G4cerr << "G4PersistentGeomMan: Error in storing geometry." << G4endl;
    return false;
  }

}

G4bool G4PersistentGeomMan::Retrieve( HepDbApplication* dbApp,
                                      G4VPhysicalVolume*& theWorld)
{
  theWorld = 0;

  G4VMaterialMap* f_Materialmap = G4VMaterialMap::GetMaterialMap();
  if(f_verboseLevel>0)
  {
    if(!f_MaterialMap)
      G4cout << "G4PersistentGeomMan::Retrieve"
             << " -- Material Map not provided by user code." << G4endl;
  }

  // ----------------------- Objectivity Feature ----------------------- //
  // Load the geometry object map from the database with read-only mode
  if ( ! f_GeomMap.lookupObj( f_DB, f_GeomMapScopeName, oocRead ) )
  {
    G4cerr << "G4PersistentGeomMan::Retrieve"
           << " -- error in looking up the geometry object map: "
           << f_GeomMapScopeName << G4endl;
    return false;
  }
  // ----------------------- Objectivity Feature ----------------------- //

  if ( f_GeomMap == 0 )
  {
    G4cerr << "G4PersistentGeomMan::Retrieve()"
           << " -- error in searching the geometry object map: "
           << f_GeomMapScopeName << G4endl;
    return false;
  }

  if(f_verboseLevel>2)
    G4cout << "G4PersistentGeomMan::Retrieve()"
           << " -- f_GeomMap is loaded." << G4endl;

  // initialize the transient part of the geometry object map
  f_GeomMap->InitTransientMap();

  if(f_verboseLevel>2)
    G4cout << "G4PersistentGeomMan::Retrieve(G4VPhysicalVolume*)"
           << " -- f_GeomMap is initialized." << G4endl;

  // loop through the number of solids in the database
  G4int nSolids = f_GeomMap->GetNoSolids();
  for(G4int i=0;i<nSolids;i++)
  {
    // load persistent solid object from DB
    HepRef(G4PVSolid) persSolid = f_GeomMap->GetPSolid(i);

    // make transient version of G4VSolid from the persistent solid
    G4VSolid* theSolid = MakeTransientObject(persSolid);
  }

  if(f_verboseLevel>2)
    G4cout << "G4PersistentGeomMan::Retrieve(G4VPhysicalVolume*)"
           << " -- Transient Solid is created for nSolids:" 
           << nSolids << G4endl;

  // Locate the persistent world physics volume
  HepRef(G4PVPhysicalVolume) persWorld = f_GeomMap->GetWorldVolume();

  if ( persWorld == 0 )
  {
    G4cerr << "G4PersistentGeomMan::Retrieve"
           << " -- geometry object map does not have world volume."
           << G4endl;
    return false;
  }

  // start the recursive data member copy
  // to create the transient geometry object tree

  HepRef(G4PVPhysicalVolume) persMother = 0;
  theWorld = MakeTransientObject( persWorld, persMother );
  
  if ( theWorld == 0 )
  {
    G4cerr << "G4PersistentGeomMan::Retrieve"
           << " -- transient world volume is empty."
           << G4endl;
    return false;
  }

  return true;
}

//----------------------------------------------------------------------------
HepRef(G4PVPhysicalVolume) G4PersistentGeomMan::MakePersistentObject (
                                 G4VPhysicalVolume* thePhysVol )
{
  // check the depth of the recursive call to this method
  if ( ++f_nRecursive > G4_PHYS_VOLUME_DEPTH_MAX )
  {
     G4cerr << "G4PersistentGeomMan::MakePersistentObject" <<
      " -- Physical Volume recursive depth reached to G4_PHYS_VOLUME_DEPTH_MAX."
       << G4endl;
    return 0;
  }

  if( f_verboseLevel>2 && f_nRecursive % 1000 == 0 )
    G4cout << "G4PersistentGeomMan: "
           << "Physical Volume recursive depth is "
           << f_nRecursive << G4endl;

  // check if the persistent version of thePhysVol exists
  HepRef(G4PVPhysicalVolume) persPhysVol = f_GeomMap->LookUp(thePhysVol);
  if (  persPhysVol == 0 )
  {
    // Get the Logical Volume in thePhysVol
    G4LogicalVolume* theLogVol = thePhysVol->GetLogicalVolume();

    // make the persistent logical volume object from "theLogVol"
    HepRef(G4PLogicalVolume) persLogVol =
                               MakePersistentObject( theLogVol );
#ifdef G4DEBUG
    assert( persLogVol != 0 );
#endif

    // Construct the persistent version of thePhysVol for each type of volume
    switch (VolumeType(thePhysVol))
      {
      case kNormal:
        persPhysVol = new(f_container)
                          G4PPVPlacement( thePhysVol, persLogVol);
        break;
      case kReplica:
        persPhysVol = new(f_container)
                          G4PPVReplica( thePhysVol, persLogVol);
        break;
      case kParameterised:
        persPhysVol = new(f_container)
                          G4PPVReplica( thePhysVol, persLogVol);
//                          G4PPVParameterized( thePhysVol, persLogVol);
        break;
      }
#ifdef G4DEBUG
    assert( persPhysVol != 0 );
#endif

    // register persPhysVol to the geometry object lookup table
    f_GeomMap->Add( thePhysVol, persPhysVol );

  }  // end of if f_GeomMap->LookUp(thePhysVol)

  return persPhysVol;
}

HepRef(G4PLogicalVolume) G4PersistentGeomMan::MakePersistentObject (
                                             G4LogicalVolume* theLogVol )

{
  // check if the persistent version of theLogVol exists
  HepRef(G4PLogicalVolume) persLogVol = f_GeomMap->LookUp(theLogVol);
  if ( persLogVol == 0 )
  {
    // get the solid of the transient logical volume
    G4VSolid* theSolid = theLogVol->GetSolid();

    // make the persistent solid object from "theSolid"
    HepRef(G4PVSolid) persSolid = MakePersistentObject(theSolid);   

    // Construct the persistent version of theLogVol
    persLogVol = new(f_container) G4PLogicalVolume(theLogVol, persSolid);
#ifdef G4DEBUG
    assert( persLogVol != 0 );
#endif

    // register persLogVol to the geometry object lookup table
    f_GeomMap->Add( theLogVol, persLogVol );

    // loop through the  daughter physical volumes in theLogVol
    G4int nDaughters = theLogVol->GetNoDaughters();
    for(G4int i=0;i<nDaughters;i++)
    {
      // Get the i-th daughter physics volume
      G4VPhysicalVolume* dPhysVol = theLogVol->GetDaughter(i);

      // make the persistent physics volume object from "dPhysVol"
      HepRef(G4PVPhysicalVolume) persDaughterPhysVol =
            MakePersistentObject( dPhysVol );

      // add persistent daughter volume to the persistent logical volume
      if ( persDaughterPhysVol != 0 )
        persLogVol->AddDaughter( persDaughterPhysVol );

    } // end of for(G4int i=0;i<nDaughters;i++)

  } // end of if f_GeomMap->LookUp(theLogVol)

  return persLogVol;
}

HepRef(G4PVSolid) G4PersistentGeomMan::MakePersistentObject (
                          G4VSolid* theSolid )
{
  // check if the persistent version of theSolid exists
  HepRef(G4PVSolid) persSolid = f_GeomMap->LookUp(theSolid);
  if ( persSolid == 0 )
  {
    G4GeometryType theSolidType = theSolid->GetEntityType();

    if(f_verboseLevel>3)
      G4cout << "G4PersistentGeomMan::MakePersistentObject (G4VSolid) "
             << "-- theSolidType is " << theSolidType << G4endl;

    // Construct persistent and concrete solid object
    // according to the entity type of theSolid

    // Pure CSG
    if            ( theSolidType == "G4Box")
      persSolid  = new(f_container) G4PBox  ( (G4Box*)theSolid );
    else if      ( theSolidType == "G4Cons")
      persSolid = new(f_container) G4PCons  ( (G4Cons*)theSolid );
    else if      ( theSolidType == "G4Para")
      persSolid = new(f_container) G4PPara  ( (G4Para*)theSolid );
    else if      ( theSolidType == "G4Sphere")
      persSolid = new(f_container) G4PSphere( (G4Sphere*)theSolid );
    else if      ( theSolidType == "G4Torus")
      persSolid = new(f_container) G4PTorus ( (G4Torus*)theSolid );
    else if      ( theSolidType == "G4Trap")
      persSolid = new(f_container) G4PTrap  ( (G4Trap*)theSolid );
    else if      ( theSolidType == "G4Trd")
      persSolid = new(f_container) G4PTrd   ( (G4Trd*)theSolid );
    else if      ( theSolidType == "G4Tubs")
      persSolid = new(f_container) G4PTubs  ( (G4Tubs*)theSolid );

   // Specific CSG
   //  else if      ( theSolidType == "G4EllipticalTube")
   //    persSolid = new(f_container) G4PEllipticalTube
   //                              ( (G4EllipticalTube*)theSolid );
    else if      ( theSolidType == "G4Hype")
      persSolid = new(f_container) G4PHype
                                 ( (G4Hype*)theSolid );
   //  else if      ( theSolidType == "G4Polycone")
   //    persSolid = new(f_container) G4PPolycone
   //                               ( (G4Polycone*)theSolid );
   //  else if      ( theSolidType == "G4Polyhedra")
   //    persSolid = new(f_container) G4PPolyhedra
   //                               ( (G4Polyhedra*)theSolid );

    // Boolian Solids
    else if ( theSolidType == "G4UnionSolid"        || 
              theSolidType == "G4SubtractionSolid"  ||
              theSolidType == "G4IntersectionSolid" )
    {
      // This is Boolean solid, so get the constituent solids
      G4BooleanSolid* boolSolid = (G4BooleanSolid*) theSolid;
      G4VSolid* aSolidA = boolSolid->GetConstituentSolid(0);
      G4VSolid* aSolidB = boolSolid->GetConstituentSolid(1);
#ifdef G4DEBUG
      assert(aSolidA!=0);
      assert(aSolidB!=0);
#endif
      HepRef(G4PVSolid) persSolidA = MakePersistentObject(aSolidA);
      HepRef(G4PVSolid) persSolidB = MakePersistentObject(aSolidB);

      const G4String solidName = theSolid->GetName();

      // Create persistent boolean solid according to its type
      if ( theSolidType == "G4UnionSolid")
      {
        persSolid = new(f_container)
          G4PUnionSolid( solidName, persSolidA, persSolidB);
      }
      else if ( theSolidType == "G4SubtractionSolid")
      {
        persSolid = new(f_container) 
           G4PSubtractionSolid( solidName, persSolidA, persSolidB);
      }
      else if ( theSolidType == "G4IntersectionSolid")
      {
        persSolid = new(f_container)
           G4PIntersectionSolid( solidName, persSolidA, persSolidB);
      }

    }
    else if ( theSolidType == "G4DisplacedSolid")
    {
      // This is a displaced solid, so get the constituent solid and transform
      G4DisplacedSolid* dispSolid = (G4DisplacedSolid*) theSolid;
      G4VSolid* movedSolid = dispSolid->GetConstituentMovedSolid();
#ifdef G4DEBUG
      assert(movedSolid!=0);
#endif
      HepRef(G4PVSolid) persMovedSolid = MakePersistentObject(movedSolid);

      HepRef(G4PAffineTransform) pTransform = new(f_container)
             G4PAffineTransform( dispSolid->GetDirectTransform() );

      persSolid = new(f_container)
         G4PDisplacedSolid( persMovedSolid, pTransform );
    }
#ifdef G4DEBUG
    else
    {
      G4cerr << " Warning (G4PersistentGeomMan): The solid type ("
             << theSolidType << ") is not supported yet." << G4endl;
    }
#endif

    // register persSolid to the geometry object lookup table
    if ( persSolid != 0 )
      f_GeomMap->Add( theSolid, persSolid );

  } // end of if f_GeomMap->LookUp(theSolid)

  return persSolid;
}

//----------------------------------------------------------------------------

G4VPhysicalVolume* G4PersistentGeomMan::MakeTransientObject (
                          HepRef(G4PVPhysicalVolume) persPhysVol,
                          HepRef(G4PVPhysicalVolume) persMotherVol )
{
  // check the depth of the recursive call to this method
  if ( ++f_tRecursive > G4_PHYS_VOLUME_DEPTH_MAX )
  {
     G4cerr << "G4PersistentGeomMan::MakeTransientObject" <<
      " -- Physical Volume recursive depth reached to G4_PHYS_VOLUME_DEPTH_MAX."
       << G4endl;
    return 0;
  }

  if( f_verboseLevel>2 && f_tRecursive % 1000 == 0 )
    G4cout << "G4PersistentGeomMan: "
           << "Physical Volume recursive depth is "
           << f_tRecursive << G4endl;

  // check if the transient version of persPhysVol exists
  G4VPhysicalVolume* thePhysVol = f_GeomMap->LookUp(persPhysVol);
  if (  thePhysVol == 0 )
  {
    // Get the Logical Volume in persPhysVol
    HepRef(G4PLogicalVolume) persLogVol = persPhysVol->GetLogicalVolume();

    // make the transient logical volume object from "persPhysVol"
    G4LogicalVolume* theLogVol =
         MakeTransientObject( persLogVol, persPhysVol );
#ifdef G4DEBUG
    assert( theLogVol != 0 );
#endif

    // Construct the transient version of persPhysVol
    G4VPhysicalVolume*
        thePhysVol = persPhysVol->MakeTransientObject(
                         theLogVol, f_GeomMap->LookUp(persMotherVol) );
#ifdef G4DEBUG
    assert( thePhysVol != 0 );
#endif

    // register thePhysVol to the geometry object lookup table
    f_GeomMap->Add( persPhysVol, thePhysVol );

  }  // end of if f_GeomMap->LookUp(persPhysVol)

  return thePhysVol;
}

G4LogicalVolume* G4PersistentGeomMan::MakeTransientObject (
                          HepRef(G4PLogicalVolume) persLogVol,
                          HepRef(G4PVPhysicalVolume) persMotherVol )
{
  // check if the transient version of persLogVol exists
  G4LogicalVolume* theLogVol = f_GeomMap->LookUp(persLogVol);
  if ( theLogVol == 0 )
  {
    // get the solid of the persistent logical volume
    HepRef(G4PVSolid) persSolid = persLogVol->GetSolid();
#ifdef G4DEBUG
    assert( persSolid != 0 );
#endif

    // check if the transient version of persSolid exists
    G4VSolid* theSolid = f_GeomMap->LookUp(persSolid);

#ifdef G4DEBUG
    if ( theSolid == 0 )
    {
      G4cerr << "G4PersistentGeomMan::MakeTransientObject" <<
       " -- transient Solid not found for the persistent Solid" << G4endl;
    }
#endif

    // Lookup the material from name
    G4Material* theMaterial = 0;
    if( f_MaterialMap )
    {
      theMaterial = f_MaterialMap->LookUp(persLogVol->GetMaterialName());
    }

    // Construct the persistent version of theLogVol with theSolid
    theLogVol = persLogVol->MakeTransientObject(theSolid, theMaterial);
#ifdef G4DEBUG
    assert( theLogVol != 0 );
#endif

    // register theLogVol to the geometry object lookup table
    f_GeomMap->Add( persLogVol, theLogVol );

    // loop through the  daughter physical volumes in persLogVol
    G4int nDaughters = persLogVol->GetNoDaughters();
    for(G4int i=0;i<nDaughters;i++)
    {
      // Get the i-th daughter physics volume
      HepRef(G4PVPhysicalVolume) persDaughterPhysVol
                        = persLogVol->GetDaughter(i);

      // make the transient physics volume object from "persDaughterPhysVol"
      G4VPhysicalVolume* theDaughterPhysVol =
            MakeTransientObject( persDaughterPhysVol, persMotherVol );

      // add transient daughter volume to the transient logical volume
      if ( theDaughterPhysVol != 0 )
        theLogVol->AddDaughter( theDaughterPhysVol );

    } // end of for(G4int i=0;i<nDaughters;i++)

  } // end of if f_GeomMap->LookUp(persLogVol)

  return theLogVol;
}

G4VSolid* G4PersistentGeomMan::MakeTransientObject (
                          HepRef(G4PVSolid) persSolid )
{
  // check if the persistent version of theSolid exists
  G4VSolid* transSolid = f_GeomMap->LookUp(persSolid);
  if ( transSolid == 0 )
  {
    G4GeometryType theSolidType = persSolid->GetEntityType();

    // Construct transient and concrete solid object
    // according to the entity type of theSolid
    if ( theSolidType == "G4UnionSolid"        || 
         theSolidType == "G4SubtractionSolid"  ||
         theSolidType == "G4IntersectionSolid" )
    {
      // This is Boolean solid, so get the constituent solids
      HepRef(G4PBooleanSolid) boolSolid = (HepRef(G4PBooleanSolid)) persSolid;
      HepRef(G4PVSolid) aSolidA = boolSolid->GetConstituentSolid(0);
      HepRef(G4PVSolid) aSolidB = boolSolid->GetConstituentSolid(1);
#ifdef G4DEBUG
      assert(aSolidA!=0);
      assert(aSolidB!=0);
#endif
      G4VSolid* transSolidA = MakeTransientObject(aSolidA);
      G4VSolid* transSolidB = MakeTransientObject(aSolidB);
      transSolid = boolSolid->MakeTransientBooleanSolid
                                (transSolidA, transSolidB);
    }
    else if ( theSolidType == "G4DisplacedSolid")
    {
      // This is a displaced solid, so get the constituent solid and 
      // its transform
      HepRef(G4PDisplacedSolid) dispSolid =
                            (HepRef(G4PDisplacedSolid)) persSolid;
      HepRef(G4PVSolid) movedSolid = dispSolid->GetConstituentMovedSolid();
#ifdef G4DEBUG
      assert(movedSolid!=0);
#endif
      G4VSolid* transMovedSolid = MakeTransientObject(movedSolid);
      transSolid = dispSolid->MakeTransientDisplacedSolid( transMovedSolid );
    }
    else
    {
      // This is a simple solid
      transSolid = persSolid->MakeTransientObject();
    }
#ifdef G4DEBUG
    assert( transSolid != 0 );
#endif
    f_GeomMap->Add( persSolid, transSolid );
  }

  return transSolid;
}

//----------------------------------------------------------------------------

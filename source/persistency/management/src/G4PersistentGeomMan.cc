// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PersistentGeomMan.cc,v 1.2 1999/03/29 16:05:27 morita Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
// class G4PersistentGeomMan 
//
// Implementation for concrete G4PersistentGeomMan.
//
// History:
// 98.10.30 Y.Morita  Splited from G4PersistencyManager

#include "G4PersistentGeomMan.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4PVSolid.hh"
#include "G4PBox.hh"
#include "G4PCons.hh"
#include "G4PHype.hh"
#include "G4PPara.hh"
#include "G4PSphere.hh"
#include "G4PTorus.hh"
#include "G4PTrap.hh"
#include "G4PTrd.hh"
#include "G4PTubs.hh"

#include "G4PLogicalVolume.hh"
#include "G4PPVPlacement.hh"
#include "G4PPVReplica.hh"
#include "G4PPVParameterised.hh"

#include "G4ios.hh"

G4PersistentGeomMan::G4PersistentGeomMan()
{
  f_GeomMapScopeName = "Geant4 Geometry Object Map";
}

G4PersistentGeomMan::~G4PersistentGeomMan()
{
}

//----------------------------------------------------------------------------

G4bool G4PersistentGeomMan::Store( HepDbApplication* dbApp,
                                   const G4VPhysicalVolume* theWorld)
{
  void* aWorld = (void*) theWorld;
  G4bool fStatus = false;
  const G4String theGeometryName = "G4EXAMPLE test geometry";

  // clear the recursive call counter
  f_nRecursive = 0;
  f_tRecursive = 0;

  // Start Updating the geometry database
  dbApp->startUpdate();

  // Open or create Geometry DB
  const HepDatabaseRef f_GeomDB = dbApp->db("Geometry");

  if ( f_GeomDB == NULL )
  {
     G4cerr << "G4PersistentGeomMan::Store" <<
      " -- error in opening or creating the geometry database." << endl;
  }

  // create a new container for geometry collection in this database
  f_GeomContainer = dbApp->container("GeometryContainer"); 

  if ( f_GeomContainer == NULL )
  {
    G4cerr << "could not find or create f_GeomContainer in the database" << endl;
  }

  // ----------------------- Objectivity Feature ----------------------- //
  // Look for the geometry object map in the database and create it
  // if it does not exist
  if ( ! f_GeomMap.lookupObj( f_GeomDB, f_GeomMapScopeName, oocUpdate ) )
  {

    // Prepare lookup table for the geometry object tree
    f_GeomMap =
       new(f_GeomContainer) G4PGeometryObjectMap( theGeometryName );

    // Add ScopeName to f_GeomMap for future reference (*Objy feature)
    if ( ! f_GeomMap.nameObj( f_GeomDB, f_GeomMapScopeName ) )
    {
       G4cerr << "G4PersistentGeomMan::Store" <<
        " -- error in naming the geometry object map: " <<
       f_GeomMapScopeName << endl;
    }

  }
  else
  {
     G4cerr << "G4PersistentGeomMan::Store" <<
      " -- the object map already defined in the geometry database." << endl;
     G4cerr << "  Geometry not stored." << endl;
     dbApp->abort();
     return false;
  }
  // ----------------------- Objectivity Feature ----------------------- //

  // make the persistent world volume object from "theWorld"
  HepRef(G4PVPhysicalVolume) persWorld =
                  MakePersistentObject( (G4VPhysicalVolume*)aWorld );

  if ( persWorld != NULL )
  {
    // set the persistent World Volume to f_GeomMap
    f_GeomMap->SetWorldVolume( persWorld );

    // commit the changes to the geometry database
    dbApp->commit();
    fStatus = true;
  }
  else
  {
    dbApp->abort();
    fStatus = false;
  }

  return fStatus;
}

G4bool G4PersistentGeomMan::Retrieve( HepDbApplication* dbApp,
                                      G4VPhysicalVolume*& theWorld)
{
  const G4PString theGeometryName = "G4EXAMPLE test geometry";
  G4bool fStatus = false;

  theWorld = NULL;

  // Start "Read" transactions to the database
  dbApp->startRead();

  // Open existing Geometry database
  f_GeomDB = dbApp->db("Geometry");

  if ( f_GeomDB == NULL )
  {
     G4cerr << "G4PersistentGeomMan::Retrieve" <<
      " -- error in searching the geometry database" << endl;
  }

  // ----------------------- Objectivity Feature ----------------------- //
  // Load the geometry object map from the database with read-only mode
  if ( ! f_GeomMap.lookupObj( f_GeomDB, f_GeomMapScopeName, oocRead ) )
  {
     G4cerr << "G4PersistentGeomMan::Retrieve" <<
      " -- error in looking up the geometry object map: " <<
     f_GeomMapScopeName << endl;
  }
  // ----------------------- Objectivity Feature ----------------------- //

  if ( f_GeomMap == NULL )
  {
     G4cerr << "G4PersistentGeomMan::Retrieve" <<
      " -- error in searching the geometry object map: " <<
     f_GeomMapScopeName << endl;
  }

#ifdef G4PERSISTENCY_DEBUG
  G4cout << "G4PersistentGeomMan::Retrieve(G4VPhysicalVolume*)"
       << " -- f_GeomMap is loaded." << endl;
#endif

  // initialize the transient part of the geometry object map
  f_GeomMap->InitTransientMap();

#ifdef G4PERSISTENCY_DEBUG
  G4cout << "G4PersistentGeomMan::Retrieve(G4VPhysicalVolume*)"
       << " -- f_GeomMap is initialized." << endl;
#endif

  // loop through the number of solids in the database
  G4int nSolids = f_GeomMap->GetNoSolids();
  for(G4int i=0;i<nSolids;i++)
  {
    // load persistent solid object from DB
    HepRef(G4PVSolid) persSolid = f_GeomMap->GetPSolid(i);

    // make transient version of G4VSolid from the persistent solid
    G4VSolid* theSolid = MakeTransientObject(persSolid);
  }

#ifdef G4PERSISTENCY_DEBUG
  G4cout << "G4PersistentGeomMan::Retrieve(G4VPhysicalVolume*)"
       << " -- Transient Solid is created for nSolids:" 
       << nSolids << endl;
#endif

  // Locate the persistent world physics volume
  HepRef(G4PVPhysicalVolume) persWorld = f_GeomMap->GetWorldVolume();

  if ( persWorld == NULL )
  {
     G4cerr << "G4PersistentGeomMan::Retrieve" <<
      " -- geometry object map does not have world volume." << endl;
  }
  else
  {
    // start the recursive data member copy
    // to create the transient geometry object tree
    HepRef(G4PVPhysicalVolume) persMother = NULL;
    G4VPhysicalVolume* theWorld =
                       MakeTransientObject( persWorld, persMother );
  
    if ( theWorld != NULL )
      fStatus = true;
  }
  // Close the transaction on geometry retrive
  dbApp->commit();

  return fStatus;
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
       << endl;
    return NULL;
  }
#ifdef G4PERSISTENCY_DEBUG
  if ( f_nRecursive % 1000 == 0 )
  {
     G4cout << "G4PersistentGeomMan::MakePersistentObject" <<
       " -- Physical Volume recursive depth is " << f_nRecursive << endl;
  }
#endif // G4PERSISTENCY_DEBUG

  // check if the persistent version of thePhysVol exists
  HepRef(G4PVPhysicalVolume) persPhysVol = f_GeomMap->LookUp(thePhysVol);
  if (  persPhysVol == NULL )
  {
    // Get the Logical Volume in thePhysVol
    G4LogicalVolume* theLogVol = thePhysVol->GetLogicalVolume();

    // make the persistent logical volume object from "theLogVol"
    HepRef(G4PLogicalVolume) persLogVol =
                               MakePersistentObject( theLogVol );
    assert( persLogVol != NULL );

    // Construct the persistent version of thePhysVol for each type of volume
    switch (VolumeType(thePhysVol))
      {
      case kNormal:
        persPhysVol = new(f_GeomContainer)
                          G4PPVPlacement( thePhysVol, persLogVol);
        break;
      case kReplica:
        persPhysVol = new(f_GeomContainer)
                          G4PPVReplica( thePhysVol, persLogVol);
        break;
      case kParameterised:
        persPhysVol = new(f_GeomContainer)
                          G4PPVReplica( thePhysVol, persLogVol);
//                          G4PPVParameterized( thePhysVol, persLogVol);
        break;
      }
    assert( persPhysVol != NULL );

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
  if ( persLogVol == NULL )
  {
    // get the solid of the transient logical volume
    G4VSolid* theSolid = theLogVol->GetSolid();

    // make the persistent solid object from "theSolid"
    HepRef(G4PVSolid) persSolid = MakePersistentObject(theSolid);   

    // Construct the persistent version of theLogVol
    persLogVol = new(f_GeomContainer) G4PLogicalVolume(theLogVol, persSolid);
    assert( persLogVol != NULL );

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
      if ( persDaughterPhysVol != NULL )
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
  if ( persSolid == NULL )
  {
    G4GeometryType theSolidType = theSolid->GetEntityType();

#ifdef G4PERSISTENCY_DEBUG
  G4cout << "G4PersistentGeomMan::MakePersistentObject (G4VSolid) "
       << "-- theSolidType is " << theSolidType << endl;
#endif

    // Construct persistent and concrete solid object
    // according to the entity type of theSolid
    if      ( theSolidType == "G4Box")
      persSolid  = new(f_GeomContainer) G4PBox( (G4Box*)theSolid );
    else if ( theSolidType == "G4Cons")
      persSolid = new(f_GeomContainer) G4PCons( (G4Cons*)theSolid );
    else if ( theSolidType == "G4Hype")
      persSolid = new(f_GeomContainer) G4PHype( (G4Hype*)theSolid );
    else if ( theSolidType == "G4Para")
      persSolid = new(f_GeomContainer) G4PPara( (G4Para*)theSolid );
    else if ( theSolidType == "G4Sphere")
      persSolid = new(f_GeomContainer) G4PSphere( (G4Sphere*)theSolid );
    else if ( theSolidType == "G4Torus")
      persSolid = new(f_GeomContainer) G4PTorus( (G4Torus*)theSolid );
    else if ( theSolidType == "G4Trap")
      persSolid = new(f_GeomContainer) G4PTrap( (G4Trap*)theSolid );
    else if ( theSolidType == "G4Trd")
      persSolid = new(f_GeomContainer) G4PTrd( (G4Trd*)theSolid );
    else if ( theSolidType == "G4Tubs")
      persSolid = new(f_GeomContainer) G4PTubs( (G4Tubs*)theSolid );

    // register persSolid to the geometry object lookup table
    if ( persSolid != NULL )
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
       << endl;
    return NULL;
  }
#ifdef G4PERSISTENCY_DEBUG
  if ( f_tRecursive % 1000 == 0 )
  {
     G4cout << "G4PersistentGeomMan::MakeTransientObject" <<
       " -- Physical Volume recursive depth is " << f_tRecursive << endl;
  }
#endif // G4PERSISTENCY_DEBUG

  // check if the transient version of persPhysVol exists
  G4VPhysicalVolume* thePhysVol = f_GeomMap->LookUp(persPhysVol);
  if (  thePhysVol == NULL )
  {
    // Get the Logical Volume in persPhysVol
    HepRef(G4PLogicalVolume) persLogVol = persPhysVol->GetLogicalVolume();

    // make the transient logical volume object from "persPhysVol"
    G4LogicalVolume* theLogVol =
         MakeTransientObject( persLogVol, persPhysVol );
    assert( theLogVol != NULL );

    // Construct the transient version of persPhysVol
    G4VPhysicalVolume*
        thePhysVol = persPhysVol->MakeTransientObject(
                         theLogVol, f_GeomMap->LookUp(persMotherVol) );
    assert( thePhysVol != NULL );

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
  if ( theLogVol == NULL )
  {
    // get the solid of the persistent logical volume
    HepRef(G4PVSolid) persSolid = persLogVol->GetSolid();
    assert( persSolid != NULL );

    // check if the transient version of persSolid exists
    G4VSolid* theSolid = f_GeomMap->LookUp(persSolid);

    if ( theSolid == NULL )
    {
      G4cerr << "G4PersistentGeomMan::MakeTransientObject" <<
       " -- transient Solid not found for the persistent Solid" << endl;
    }

    // Construct the persistent version of theLogVol with theSolid
    theLogVol = persLogVol->MakeTransientObject(theSolid);
    assert( theLogVol != NULL );

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
      if ( theDaughterPhysVol != NULL )
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
  if ( transSolid == NULL )
  {
#ifdef G4PERSISTENCY_DEBUG
    G4GeometryType theSolidType = persSolid->GetEntityType();
    G4cout << "G4PersistentGeomMan::MakePersistentObject (G4PVSolid) "
           << "-- theSolidType is " << theSolidType << endl;
#endif
    transSolid = persSolid->MakeTransientObject();
    assert( transSolid != NULL );
    f_GeomMap->Add( persSolid, transSolid );
  }

  return transSolid;
}

//----------------------------------------------------------------------------

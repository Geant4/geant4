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
// $Id: G4PGeometryObjectMap.cc,v 1.7.2.1 2001/06/28 19:11:28 gunter Exp $
// GEANT4 tag $Name:  $
//
// class G4PGeometryObjectMap 
//
// Implementation of a class responsible for keeping track of
// the geometry object map.
//
// History:
// 98.06.20 Y.Morita  Initial version

#include <assert.h>
#include "G4PGeometryObjectMap.hh"

// forward declarations
#include "G4PVPhysicalVolume.hh"
#include "G4PLogicalVolume.hh"
#include "G4PVSolid.hh"

#include "G4VPhysVolRefArray.hh"
#include "G4LogVolRefArray.hh"
#include "G4VSolidRefArray.hh"

G4PGeometryObjectMap::G4PGeometryObjectMap()
 : noPhysVol(0),noLogVol(0),noSolids(0)
{
  transPhysVolPtrs = new G4VPhysVolRefArray;
  transLogVolPtrs  = new G4LogVolRefArray;
  transSolidPtrs   = new G4VSolidRefArray;
}

G4PGeometryObjectMap::G4PGeometryObjectMap( const G4String theGeometryName )
 : noPhysVol(0),noLogVol(0),noSolids(0)
{
  transPhysVolPtrs = new G4VPhysVolRefArray;
  transLogVolPtrs  = new G4LogVolRefArray;
  transSolidPtrs   = new G4VSolidRefArray;
}

G4PGeometryObjectMap::~G4PGeometryObjectMap()
{
  delete transPhysVolPtrs;
  delete transLogVolPtrs;
  delete transSolidPtrs;
}

// -------- LookUp() and Add() functions for each class -------- //

HepRef(G4PVPhysicalVolume) G4PGeometryObjectMap::LookUp(
                                      G4VPhysicalVolume* inGeomObj )
{
  assert(inGeomObj != 0);
  for(G4int i=0;i<noPhysVol;i++)
  {
    if( transPhysVolPtrs->Get(i) == inGeomObj )
    {
      return persPhysVolPtrs[i];
    }
  }

  return 0;
}

void G4PGeometryObjectMap::Add(         G4VPhysicalVolume*  inGeomObj,
                                HepRef(G4PVPhysicalVolume) outGeomObj )
{
  assert(inGeomObj != 0);
  HepRef(G4PVPhysicalVolume) aGeomObj = LookUp( inGeomObj );
  if( aGeomObj == 0 )
  {
//    assert( inGeomObj == aGeomObj );
    noPhysVol++;
    transPhysVolPtrs->Resize(noPhysVol);
    transPhysVolPtrs->Insert(noPhysVol-1, inGeomObj);
    persPhysVolPtrs.insert_element(outGeomObj);
  }
}

// ------------------------------------------------------------- //

G4VPhysicalVolume* G4PGeometryObjectMap::LookUp(
                               HepRef(G4PVPhysicalVolume) inGeomObj )
{
  assert(inGeomObj != 0);
  for(G4int i=0;i<noPhysVol;i++)
  {
    if( persPhysVolPtrs[i] == inGeomObj )
    {
      return transPhysVolPtrs->Get(i);
    }
  }

  return 0;
}

void G4PGeometryObjectMap::Add( HepRef(G4PVPhysicalVolume)  inGeomObj,
                                        G4VPhysicalVolume* outGeomObj )
{
  assert(inGeomObj != 0);
  for(G4int i=0;i<noPhysVol;i++)
  {
    if( persPhysVolPtrs[i] == inGeomObj )
    {
      transPhysVolPtrs->Insert(i, outGeomObj);
      break;
    }
    G4cerr << "G4PGeometryObjectMap::Add -- transient solid not assigned" << G4endl;
  }
}

// ------------------------------------------------------------- //

HepRef(G4PLogicalVolume) G4PGeometryObjectMap::LookUp(
                                      G4LogicalVolume* inGeomObj )
{
  assert(inGeomObj != 0);
  for(G4int i=0;i<noLogVol;i++)
  {
    if( transLogVolPtrs->Get(i) == inGeomObj )
    {
      return persLogVolPtrs[i];
    }
  }

  return 0;
}

void G4PGeometryObjectMap::Add(         G4LogicalVolume*  inGeomObj,
                                HepRef(G4PLogicalVolume) outGeomObj )
{
  assert(inGeomObj != 0);
  HepRef(G4PLogicalVolume) aGeomObj = LookUp( inGeomObj );
  if( aGeomObj == 0 )
  {
//    assert( inGeomObj == aGeomObj );
    noLogVol++;
    transLogVolPtrs->Resize(noLogVol);
    transLogVolPtrs->Insert(noLogVol-1, inGeomObj);
    persLogVolPtrs.insert_element(outGeomObj);
  }
}

// ------------------------------------------------------------- //

G4LogicalVolume* G4PGeometryObjectMap::LookUp(
                               HepRef(G4PLogicalVolume) inGeomObj )
{
  assert(inGeomObj != 0);
  for(G4int i=0;i<noLogVol;i++)
  {
    if( persLogVolPtrs[i] == inGeomObj )
    {
      return transLogVolPtrs->Get(i);
    }
  }

  return 0;
}

void G4PGeometryObjectMap::Add( HepRef(G4PLogicalVolume)  inGeomObj,
                                        G4LogicalVolume* outGeomObj )
{
  assert(inGeomObj != 0);
  for(G4int i=0;i<noLogVol;i++)
  {
    if( persLogVolPtrs[i] == inGeomObj )
    {
      transLogVolPtrs->Insert(i, outGeomObj);
      break;
    }
    G4cerr << "G4PGeometryObjectMap::Add -- transient solid not assigned" << G4endl;
  }
}

// ------------------------------------------------------------- //
HepRef(G4PVSolid) G4PGeometryObjectMap::LookUp( G4VSolid* inGeomObj )
{

#ifdef G4PERSISTENCY_DEBUG
  cout << "G4PGeometryObjectMap::LookUp(G4VSolid) -- noSolids is "
       << noSolids << G4endl;
  cout << "  G4VSolid info: " << inGeomObj << G4endl;
#endif

  assert(inGeomObj != 0);
  for(G4int i=0;i<noSolids;i++)
  {
#ifdef G4PERSISTENCY_DEBUG
    G4VSolid* tmpSolid = transSolidPtrs->Get(i);
    cout << "[" << i << "] transSolidPtrs[i]=" << tmpSolid <<G4endl;
    cout             << "   persSolidPtrs[i]=";
    HepRef(G4PVSolid) tmpPSolid = persSolidPtrs[i];
    tmpPSolid.print();
#endif

    if( transSolidPtrs->Get(i) == inGeomObj )
    {
      return persSolidPtrs[i];
    }
  }

  return 0;
}

void G4PGeometryObjectMap::Add(         G4VSolid*  inGeomObj,
                                HepRef(G4PVSolid) outGeomObj )
{
  assert(inGeomObj != 0);
  HepRef(G4PVSolid) aGeomObj = LookUp( inGeomObj );
  if( aGeomObj == 0 )
  {
//    assert( inGeomObj == aGeomObj );
    noSolids++;
    transSolidPtrs->Resize(noSolids);
    transSolidPtrs->Insert(noSolids-1, inGeomObj);
    persSolidPtrs.insert_element(outGeomObj);
  }
}

// ------------------------------------------------------------- //

G4VSolid* G4PGeometryObjectMap::LookUp(
                               HepRef(G4PVSolid) inGeomObj )
{

#ifdef G4PERSISTENCY_DEBUG
  cout << "G4PGeometryObjectMap::LookUp(G4PVSolid) -- noSolids is "
       << noSolids << G4endl;
  cout << "  G4PVSolid info:";
  inGeomObj.print(stdout);
#endif

  assert(inGeomObj != 0);
  for(G4int i=0;i<noSolids;i++)
  {
#ifdef G4PERSISTENCY_DEBUG
    cout << "[" << i << "]  persSolidPtrs[i]=";
    HepRef(G4PVSolid) tmpPSolid = persSolidPtrs[i];
    tmpPSolid.print();
    G4VSolid* tmpSolid = transSolidPtrs->Get(i);
    cout              << " transSolidPtrs[i]=" << tmpSolid <<G4endl;
#endif
    if( persSolidPtrs[i] == inGeomObj )
    {
      return transSolidPtrs->Get(i);
    }
  }

  return 0;
}

void G4PGeometryObjectMap::Add( HepRef(G4PVSolid)  inGeomObj,
                                        G4VSolid* outGeomObj )
{
  assert(inGeomObj != 0);
  for(G4int i=0;i<noSolids;i++)
  {
    if( persSolidPtrs[i] == inGeomObj )
    {
      transSolidPtrs->Insert(i, outGeomObj);
      break;
    }
    G4cerr << "G4PGeometryObjectMap::Add -- transient solid not assigned" << G4endl;
  }
}

// ------------------------------------------------------------- //

G4VSolid* G4PGeometryObjectMap::GetSolid(G4int i)
{
  if( i >= 0 && i < noSolids)
  {
    return transSolidPtrs->Get(i);
  }
  else
  {
    return 0;
  }
}

HepRef(G4PVSolid) G4PGeometryObjectMap::GetPSolid(G4int i)
{
  if( i >= 0 && i < noSolids)
  {
    return persSolidPtrs[i];
  }
  else
  {
    return 0;
  }
}

void G4PGeometryObjectMap::InitTransientMap()
{
  transPhysVolPtrs->Resize(noPhysVol);
  transLogVolPtrs->Resize(noLogVol);
  transSolidPtrs->Resize(noSolids);

#ifdef G4PERSISTENCY_DEBUG
  cout << "G4PGeometryObjectMap::InitTransientMap()" << G4endl;
  cout << " -- noPhysVol: " << noPhysVol << G4endl;
  cout << " -- transPhysVolPtrs.length: " << transPhysVolPtrs->length() << G4endl;
  cout << " -- transLogVolPtrs.length: "  << transLogVolPtrs->length()  << G4endl;
  cout << " -- transSolidPtrs.length: "   << transSolidPtrs->length()   << G4endl;
#endif

  G4int i;
  for(i=0;i<noPhysVol;i++)
  {
    transPhysVolPtrs->Insert(i, 0);
  }
  for(i=0;i<noLogVol;i++)
  {
    transLogVolPtrs->Insert(i, 0);
  }
  for(i=0;i<noSolids;i++)
  {
    transSolidPtrs->Insert(i, 0);
  }

}

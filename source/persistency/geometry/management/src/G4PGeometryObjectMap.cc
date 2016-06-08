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
// $Id: G4PGeometryObjectMap.cc,v 1.9 2001/07/11 10:02:19 gunter Exp $
// GEANT4 tag $Name: geant4-04-00 $
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
    if( (*transPhysVolPtrs)[i] == inGeomObj )
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
    transPhysVolPtrs->push_back(inGeomObj);
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
      return (*transPhysVolPtrs)[i];
    }
  }

  return 0;
}

void G4PGeometryObjectMap::Add( HepRef(G4PVPhysicalVolume)  inGeomObj,
                                        G4VPhysicalVolume* outGeomObj )
{
  assert(inGeomObj != 0);
  for(G4int i=0; i<noPhysVol; i++)
  {
    if( persPhysVolPtrs[i] == inGeomObj )
    {
      (*transPhysVolPtrs)[i]= outGeomObj;
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
    if( (*transLogVolPtrs)[i] == inGeomObj )
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
    transLogVolPtrs->push_back(inGeomObj);
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
      return (*transLogVolPtrs)[i];
    }
  }

  return 0;
}

void G4PGeometryObjectMap::Add( HepRef(G4PLogicalVolume)  inGeomObj,
                                        G4LogicalVolume* outGeomObj )
{
  assert(inGeomObj != 0);
  for(G4int i=0; i<noLogVol; i++)
  {
    if( persLogVolPtrs[i] == inGeomObj )
    {
      (*transLogVolPtrs)[i] = outGeomObj;
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
    G4VSolid* tmpSolid = (*transSolidPtrs)[i];
    cout << "[" << i << "] transSolidPtrs[i]=" << tmpSolid <<G4endl;
    cout             << "   persSolidPtrs[i]=";
    HepRef(G4PVSolid) tmpPSolid = persSolidPtrs[i];
    tmpPSolid.print();
#endif

    if( (*transSolidPtrs)[i] == inGeomObj )
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
    transSolidPtrs->push_back(inGeomObj);
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
    G4VSolid* tmpSolid = (*transSolidPtrs)[i];
    cout              << " transSolidPtrs[i]=" << tmpSolid <<G4endl;
#endif
    if( persSolidPtrs[i] == inGeomObj )
    {
      return (*transSolidPtrs)[i];
    }
  }

  return 0;
}

void G4PGeometryObjectMap::Add( HepRef(G4PVSolid)  inGeomObj,
                                        G4VSolid* outGeomObj )
{
  assert(inGeomObj != 0);
  for(G4int i=0; i<noSolids; i++)
  {
    if( persSolidPtrs[i] == inGeomObj )
    {
      (*transSolidPtrs)[i] = outGeomObj;
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
    return (*transSolidPtrs)[i];
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
  transPhysVolPtrs->resize(noPhysVol);
  transLogVolPtrs->resize(noLogVol);
  transSolidPtrs->resize(noSolids);

#ifdef G4PERSISTENCY_DEBUG
  cout << "G4PGeometryObjectMap::InitTransientMap()" << G4endl;
  cout << " -- noPhysVol: " << noPhysVol << G4endl;
  cout << " -- transPhysVolPtrs.length: " << transPhysVolPtrs->size() << G4endl;
  cout << " -- transLogVolPtrs.length: "  << transLogVolPtrs->size()  << G4endl;
  cout << " -- transSolidPtrs.length: "   << transSolidPtrs->size()   << G4endl;
#endif

  G4int i;
  for(i=0;i<noPhysVol;i++)
  {
    transPhysVolPtrs->push_back(0);
  }
  for(i=0;i<noLogVol;i++)
  {
    transLogVolPtrs->push_back(0);
  }
  for(i=0;i<noSolids;i++)
  {
    transSolidPtrs->push_back(0);
  }

}

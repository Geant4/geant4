// $Id: ExP02GeoTree.cc,v 1.1 2007-11-29 17:05:22 witoldp Exp $
// Include files

// local
#include "ExP02GeoTree.hh"

//-----------------------------------------------------------------------------
// Implementation file for class : ExP02GeoTree
//
// 2005-05-26 : Witold POKORSKI
//-----------------------------------------------------------------------------

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
ExP02GeoTree::ExP02GeoTree(  ): topV(0)
{
  eltab = G4Element::GetElementTable();
  mattab = G4Material::GetMaterialTable();
}

ExP02GeoTree::ExP02GeoTree(G4VPhysicalVolume* vol, const G4ElementTable* et, const G4MaterialTable* mt): 
  topV(vol), eltab(et), mattab(mt)
{}

//=============================================================================
// Destructor
//=============================================================================
ExP02GeoTree::~ExP02GeoTree() {} 

//=============================================================================

G4VPhysicalVolume* ExP02GeoTree::TopVol()
{
  return topV;
}

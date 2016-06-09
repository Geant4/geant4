//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: ExP02GeoTree.cc,v 1.2 2007-12-10 16:29:20 gunter Exp $
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

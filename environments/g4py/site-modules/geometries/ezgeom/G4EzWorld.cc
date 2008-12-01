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
// $Id: G4EzWorld.cc,v 1.1 2008-12-01 07:07:34 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   G4EzWorld.cc
//
//                                         2005 Q
// ====================================================================
#include "G4EzWorld.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4RunManager.hh"

G4VPhysicalVolume* G4EzWorld::world= G4EzWorld::CreateWorld();

// ====================================================================
//
// class description
//
// ====================================================================

//////////////////////
G4EzWorld::G4EzWorld()
//////////////////////
{
}

///////////////////////
G4EzWorld::~G4EzWorld()
///////////////////////
{
  world= 0;
}


//////////////////////////////////////////////////////////
G4VPhysicalVolume* G4EzWorld::CreateWorld
                   (G4double dx, G4double dy, G4double dz)
//////////////////////////////////////////////////////////
{
  // default matetial is "vacuum"
  G4Material* vacuum= G4Material::GetMaterial("_Vacuum", false);

  if(vacuum==0) {
    G4Element* elN= new G4Element("_N", "",  7.,  14.00674*g/mole);
    G4Element* elO= new G4Element("_O", "",  8.,  15.9994*g/mole);

    vacuum= new G4Material("_Vacuum", universe_mean_density, 2);
    vacuum-> AddElement(elN,  0.7);
    vacuum-> AddElement(elO,  0.3);
  }

  G4Box* sdworld= new G4Box("world", dx/2., dy/2., dz/2.);
  G4LogicalVolume* lvworld= new G4LogicalVolume(sdworld, vacuum, "word");
  G4PVPlacement* aworld= new G4PVPlacement(0, G4ThreeVector(), "world",
                                           lvworld, 0, false, 0);
  
  // vis. attributes
  G4VisAttributes* vaworld= new G4VisAttributes(G4Color(1.,1.,1.));
  vaworld-> SetForceWireframe(true);
  lvworld-> SetVisAttributes(vaworld);

  return aworld;

}


////////////////////////////////////////////////////////////
void G4EzWorld::Reset(G4double dx, G4double dy, G4double dz)
////////////////////////////////////////////////////////////
{
  delete world;
  world= CreateWorld(dx, dy, dz);

  G4RunManager* runManager= G4RunManager::GetRunManager();
  runManager-> DefineWorldVolume(world);
}


/////////////////////////////////////////////////////////////
void G4EzWorld::Resize(G4double dx, G4double dy, G4double dz)
/////////////////////////////////////////////////////////////
{
  G4Box* box= dynamic_cast<G4Box*>(world-> GetLogicalVolume()-> GetSolid());
  box-> SetXHalfLength(dx/2.);
  box-> SetYHalfLength(dy/2.);
  box-> SetZHalfLength(dz/2.);

  G4RunManager* runManager= G4RunManager::GetRunManager();
  runManager-> GeometryHasBeenModified();
}


//////////////////////////////////////////////////
void G4EzWorld::SetMaterial(G4Material* amaterial)
//////////////////////////////////////////////////
{
  world-> GetLogicalVolume()-> SetMaterial(amaterial);
}


//////////////////////////////////////////
void G4EzWorld::SetVisibility(G4bool qvis)
//////////////////////////////////////////
{
  G4VisAttributes* vaworld= const_cast<G4VisAttributes*>(
    world-> GetLogicalVolume()-> GetVisAttributes());
  vaworld-> SetVisibility(qvis);
}


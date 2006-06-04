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
// $Id: G4EzWorld.cc,v 1.2 2006-06-04 21:36:34 kmura Exp $
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
#include "G4Version.hh"

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
#if G4VERSION_NUMBER >= 800
  G4Material* vacuum= G4Material::GetMaterial("_Vacuum", false);
#else
  G4Material* vacuum= G4Material::GetMaterial("_Vacuum");
#endif

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


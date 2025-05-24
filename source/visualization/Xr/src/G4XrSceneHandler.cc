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
//
//
//
// John Allison  5th April 2001
// A template for a simplest possible graphics driver.
//?? Lines or sections marked like this require specialisation for your driver.

#include "G4XrSceneHandler.hh"

#include "G4Box.hh"
#include "G4Circle.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeModel.hh"
#include "G4Material.hh"
#include "G4Mesh.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4Polyhedron.hh"
#include "G4Polyline.hh"
#include "G4PseudoScene.hh"
#include "G4Square.hh"
#include "G4SystemOfUnits.hh"
#include "G4Text.hh"
#include "G4UnitsTable.hh"
#include "G4VNestedParameterisation.hh"
#include "G4VPhysicalVolume.hh"

// Counter for Xr scene handlers.
G4int G4XrSceneHandler::fSceneIdCount = 0;

G4XrSceneHandler::G4XrSceneHandler(G4VGraphicsSystem& system, const G4String& name)
  : G4VSceneHandler(system, fSceneIdCount++, name)
{}

void G4XrSceneHandler::AddPrimitive(const G4Polyline& polyline)
{
  std::cout << "G4XrSceneHandler::AddPrimitive(G4Polyline=" << &polyline << ")" << std::endl;
}

void G4XrSceneHandler::AddPrimitive(const G4Text& text)
{
  std::cout << "G4XrSceneHandler::AddPrimitive(G4Text=" << &text << ")" << std::endl;
}


void G4XrSceneHandler::AddPrimitive(const G4Circle& circle)
{
  std::cout << "G4XrSceneHandler::AddPrimitive(G4Circle=" << &circle << ")" << std::endl;
}

void G4XrSceneHandler::AddPrimitive(const G4Square& square)
{
  std::cout << "G4XrSceneHandler::AddPrimitive(G4Square=" << &square << ")" << std::endl;
}

void G4XrSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron)
{
  std::cout << "G4XrSceneHandler::AddPrimitive(G4Polyhedron=" << &polyhedron << ")" << std::endl;
}



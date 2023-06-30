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

#include "G4VtkQtSceneHandler.hh"

#include "G4Box.hh"
#include "G4Circle.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeModel.hh"
#include "G4Material.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4Polyhedron.hh"
#include "G4Polyline.hh"
#include "G4Square.hh"
#include "G4SystemOfUnits.hh"
#include "G4Text.hh"
#include "G4UnitsTable.hh"
#include "G4VPhysicalVolume.hh"

G4int G4VtkQtSceneHandler::fSceneIdCount = 0;
// Counter for XXX scene handlers.

G4VtkQtSceneHandler::G4VtkQtSceneHandler(G4VGraphicsSystem& system, const G4String& name)
  : G4VtkSceneHandler(system, name)
{}

void G4VtkQtSceneHandler::AddPrimitive(const G4Polyline& polyline)
{
#ifdef G4VTKDEBUG
  G4cout << fpModel->GetGlobalDescription() << " " << fpModel->GetType() << G4endl;
#endif

  // Normal base class behaviour
  G4VtkSceneHandler::AddPrimitive(polyline);
}

void G4VtkQtSceneHandler::AddPrimitive(const G4Text& text)
{
#ifdef G4VTKDEBUG
  G4cout << fpModel->GetGlobalDescription() << " " << fpModel->GetType() << G4endl;
#endif

  // Normal base class behaviour
  G4VtkSceneHandler::AddPrimitive(text);
}

void G4VtkQtSceneHandler::AddPrimitive(const G4Circle& circle)
{
#ifdef G4VTKDEBUG
  G4cout << fpModel->GetGlobalDescription() << " " << fpModel->GetType() << G4endl;
#endif

  // Normal base class behaviour
  G4VtkSceneHandler::AddPrimitive(circle);
}

void G4VtkQtSceneHandler::AddPrimitive(const G4Square& square)
{
#ifdef G4VTKDEBUG
  G4cout << fpModel->GetGlobalDescription() << " " << fpModel->GetType() << G4endl;
#endif

  // Normal base class behaviour
  G4VtkSceneHandler::AddPrimitive(square);
}

void G4VtkQtSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron)
{
#ifdef G4VTKDEBUG
  G4cout << fpModel->GetGlobalDescription() << " " << fpModel->GetType() << G4endl;
#endif

  // Normal base class behaviour
  G4VtkSceneHandler::AddPrimitive(polyhedron);
}

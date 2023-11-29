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
// Abstract interface for GEANT4 Visualization Manager.
// John Allison 19/Oct/1996.

#include "G4VVisManager.hh"

G4VVisManager::G4VVisManager ()
= default;

G4VVisManager::~G4VVisManager () = default;

G4VVisManager* G4VVisManager::fpConcreteInstance = nullptr;

G4VVisManager* G4VVisManager::GetConcreteInstance ()
{
  return fpConcreteInstance;
}

void G4VVisManager::SetConcreteInstance (G4VVisManager* man)
{
  fpConcreteInstance = man;
}

void G4VVisManager::DrawGeometry
(G4VPhysicalVolume*, const G4Transform3D&)
// Draws a geometry tree starting at the specified physical volume.
{}

void G4VVisManager::IgnoreStateChanges (G4bool)
{;}


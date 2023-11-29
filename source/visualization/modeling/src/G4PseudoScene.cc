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
// John Allison  15th May 2014
// A base class for "pseudo" scenes/graphics systems.

#include "G4PseudoScene.hh"

#include "G4Mesh.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

void G4PseudoScene::AddCompound(const G4Mesh& mesh) {
  // Catches mesh if special mesh rendering set
  ProcessVolume(*mesh.GetContainerVolume()->GetLogicalVolume()->GetSolid());
}

void G4PseudoScene::ProcessVolume (const G4VSolid& solid)
{
  G4ExceptionDescription ed;
  ed << "G4PseudoScene::ProcessVolume called for solid \"" << solid.GetName()
  << "\".\n  This is a base class - it shouldn't happen."
  << "\n  The concrete implementation has not processed this solid.";
  G4Exception("G4PseudoScene::ProcessVolume", "modeling0014",
	      FatalException, ed);
}

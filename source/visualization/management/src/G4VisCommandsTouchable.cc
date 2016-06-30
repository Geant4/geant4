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
// $Id: G4VisCommandsTouchable.cc 75794 2013-11-06 13:25:30Z allison $

// /vis/touchable commands - John Allison  14th May 2014

#include "G4VisCommandsTouchable.hh"

#include "G4UIcmdWithoutParameter.hh"
#include "G4TransportationManager.hh"
#include "G4TouchableDumpScene.hh"

G4VisCommandsTouchable::G4VisCommandsTouchable()
{
  fpCommandDump = new G4UIcmdWithoutParameter("/vis/touchable/dump",this);
  fpCommandDump->SetGuidance("Dump touchable attributes.");
  fpCommandDump->SetGuidance
  ("Use \"/vis/set/touchable\" to set current touchable.");
  fpCommandDump->SetGuidance
  ("Use \"/vis/touchable/set\" to set attributes.");
}

G4VisCommandsTouchable::~G4VisCommandsTouchable() {
  delete fpCommandDump;
}

G4String G4VisCommandsTouchable::GetCurrentValue(G4UIcommand*) {
  return "";
}

void G4VisCommandsTouchable::SetNewValue
(G4UIcommand* command,G4String)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4TransportationManager* transportationManager =
  G4TransportationManager::GetTransportationManager ();

  size_t nWorlds = transportationManager->GetNoWorlds();

  G4VPhysicalVolume* world = *(transportationManager->GetWorldsIterator());
  if (!world) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr <<
      "ERROR: G4VisCommandsTouchable::SetNewValue:"
      "\n  No world.  Maybe the geometry has not yet been defined."
      "\n  Try \"/run/initialize\""
      << G4endl;
    }
    return;
  }

  if (command == fpCommandDump)
  {
    G4bool found = false;
    std::vector<G4VPhysicalVolume*>::iterator iterWorld =
    transportationManager->GetWorldsIterator();
    for (size_t i = 0; i < nWorlds; ++i, ++iterWorld) {
      G4PhysicalVolumeModel pvModel (*iterWorld);  // Unlimited depth.
      G4ModelingParameters mp;  // Default - no culling.
      pvModel.SetModelingParameters (&mp);
      G4TouchableDumpScene dumpScene (G4cout,&pvModel,fCurrentTouchablePath);
      pvModel.DescribeYourselfTo (dumpScene);  // Initiate dump.
      if (dumpScene.IsFound()) found = true;
    }
    if (!found) {
      G4cout << "Touchable not found." << G4endl;
    }
    return;

  } else {

    if (verbosity >= G4VisManager::errors) {
      G4cerr <<
      "ERROR: G4VisCommandsTouchable::SetNewValue: unrecognised command."
      << G4endl;
    }
    return;
  }
}

// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ASCIITreeSceneHandler.cc,v 1.1 2001-04-10 15:09:33 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A scene handler to dump geometry hierarchy.
// Based on a provisional G4ASCIITreeGraphicsScene (was in modeling).

#include "G4ASCIITreeSceneHandler.hh"

#include "G4VSolid.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

G4ASCIITreeSceneHandler::G4ASCIITreeSceneHandler(G4VGraphicsSystem& system,
						 const G4String& name):
  G4VTreeSceneHandler(system, name) {}

G4ASCIITreeSceneHandler::~G4ASCIITreeSceneHandler () {}

void G4ASCIITreeSceneHandler::Dump (const G4VSolid& solid) {
  for (G4int i = 0; i < fCurrentDepth; i++ ) G4cout << "  ";
  G4cout << "\"" << fpCurrentPV -> GetName ()
	 << "\", copy no. " << fpCurrentPV -> GetCopyNo ()
	 << G4endl;
  return;
}

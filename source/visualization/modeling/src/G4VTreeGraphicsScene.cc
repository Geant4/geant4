// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VTreeGraphicsScene.cc,v 1.1 2000-05-22 07:39:41 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  18th May 2000
// An artificial graphics scene to find dump geometry hierarchy.

#include "G4VTreeGraphicsScene.hh"

#include "G4VSolid.hh"
#include "G4Vector3D.hh"
#include "G4PhysicalVolumeModel.hh"

G4VTreeGraphicsScene::G4VTreeGraphicsScene():
  fCurrentDepth                 (0),
  fpCurrentPV                   (0),
  fpCurrentLV                   (0),
  fpCurrentObjectTransformation (0)
{}

G4VTreeGraphicsScene::~G4VTreeGraphicsScene () {}

void G4VTreeGraphicsScene::EstablishSpecials
(G4PhysicalVolumeModel& pvModel) {
  pvModel.DefinePointersToWorkingSpace (&fCurrentDepth,
					&fpCurrentPV,
					&fpCurrentLV);
}

void G4VTreeGraphicsScene::Dump (const G4VSolid& solid) {
  for (G4int i = 0; i < fCurrentDepth; i++ ) G4cout << "  ";
  G4cout << "\"" << fpCurrentPV -> GetName ()
	 << "\", copy no. " << fpCurrentPV -> GetCopyNo ()
	 << G4endl;
  return;
}

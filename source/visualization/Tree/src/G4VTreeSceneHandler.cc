// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VTreeSceneHandler.cc,v 1.3 2001-06-15 07:12:37 stanaka Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A base class for a scene handler to dump geometry hierarchy.
// Based on a provisional G4VTreeGraphicsScene (was in modeling).

#include "G4VTreeSceneHandler.hh"

#include "G4VSolid.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4ModelingParameters.hh"

G4int G4VTreeSceneHandler::fSceneIdCount = 0;
// Counter for Tree scene handlers.

G4int G4VTreeSceneHandler::fSceneCount = 0;
// No. of extanct scene handlers.

G4VTreeSceneHandler::G4VTreeSceneHandler(G4VGraphicsSystem& system,
					 const G4String& name):
  G4VSceneHandler(system, fSceneIdCount++, name),
  fCurrentDepth                 (0),
  fpCurrentPV                   (0),
  fpCurrentLV                   (0),
  fpCurrentObjectTransformation (0)
{
  fSceneCount++;
}

G4VTreeSceneHandler::~G4VTreeSceneHandler () {}

void G4VTreeSceneHandler::EstablishSpecials
(G4PhysicalVolumeModel& pvModel) {
  pvModel.DefinePointersToWorkingSpace (&fCurrentDepth,
					&fpCurrentPV,
					&fpCurrentLV);
}

void G4VTreeSceneHandler::BeginModeling() {
  G4VSceneHandler::BeginModeling();  // Required: see G4VSceneHandler.hh.
  // Force culling off...
  if (fpModel) {
    fpOriginalMP = fpModel->GetModelingParameters();
    if (fpOriginalMP) {
      fpNonCullingMP = new G4ModelingParameters(*fpOriginalMP);
      fpNonCullingMP->SetCulling(false);
      ((G4VModel*)fpModel)->SetModelingParameters(fpNonCullingMP);
      // Note the deliberate casting away of const.
    }
  }
}

void G4VTreeSceneHandler::EndModeling() {
  if (fpModel && fpOriginalMP) {
    ((G4VModel*)fpModel)->SetModelingParameters(fpOriginalMP);
    delete fpNonCullingMP;
  }
  G4VSceneHandler::EndModeling();  // Required: see G4VSceneHandler.hh.
}

// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ASCIITreeSceneHandler.cc,v 1.2 2001-05-18 10:03:13 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A scene handler to dump geometry hierarchy.
// Based on a provisional G4ASCIITreeGraphicsScene (was in modeling).

#include "G4ASCIITreeSceneHandler.hh"

#include "G4ASCIITree.hh"
#include "G4VSolid.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VPVParameterisation.hh"

G4ASCIITreeSceneHandler::G4ASCIITreeSceneHandler(G4VGraphicsSystem& system,
						 const G4String& name):
  G4VTreeSceneHandler(system, name) {}

G4ASCIITreeSceneHandler::~G4ASCIITreeSceneHandler () {}

void G4ASCIITreeSceneHandler::BeginModeling () {
  fReplicaSet.clear();
}

void G4ASCIITreeSceneHandler::EndModeling () {
  fReplicaSet.clear();
}

void G4ASCIITreeSceneHandler::Dump (const G4VSolid& solid) {

  const G4ASCIITree* pSystem = (G4ASCIITree*)GetGraphicsSystem();
  const G4int verbosity =pSystem->GetVerbosity();

  if (verbosity < 1 && fReplicaSet.find(fpCurrentPV) != fReplicaSet.end()) {
    // Ignore if an already treated replica.  (Assumes model which has
    // invoked this function is a G4PhysicalVolumeModel - how can we
    // check this?)
    ((G4PhysicalVolumeModel*)fpModel)->CurtailDescent();
    return;
  }

  for (G4int i = 0; i < fCurrentDepth; i++ ) G4cout << "  ";
  G4cout << "\"" << fpCurrentPV->GetName()
	 << "\", copy no. " << fpCurrentPV->GetCopyNo();
  if (fpCurrentPV->IsReplicated()) {
    fReplicaSet.insert(fpCurrentPV);
    EAxis axis;
    G4int nReplicas;
    G4double width;
    G4double offset;
    G4bool consuming;
    fpCurrentPV->GetReplicationData(axis,nReplicas,width,offset,consuming);
    G4VPVParameterisation* pP = fpCurrentPV->GetParameterisation();
    G4cout << " (" << nReplicas;
    if (pP) {  // Parametrised volume.
      G4cout << " parametrised volumes)";
    }
    else {  // Plain replicated volume.  From geometry_guide.txt...
      G4cout << " plain replicas)";
    }
  }
  G4cout << G4endl;
  return;
}

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
//
// $Id: G4ASCIITreeSceneHandler.cc,v 1.9 2001-07-11 10:09:08 gunter Exp $
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

  G4VTreeSceneHandler::BeginModeling ();  // To re-use "culling off" code.

  const G4ASCIITree* pSystem = (G4ASCIITree*)GetGraphicsSystem();
  const G4int verbosity = pSystem->GetVerbosity();
  G4cout << "\nG4ASCIITreeSceneHandler::BeginModeling:"
    "\n  set verbosity with \"/vis/ASCIITree/verbose <verbosity>\":"
    "\n  <  10: - does not print daughters of repeated logical volumes."
    "\n         - does not repeat replicas."
    "\n  >= 10: prints all physical volumes."
    "\n  For level of detail add:"
    "\n  >=  0: prints physical volume name."
    "\n  >=  1: prints logical volume name."
    "\n  >=  2: prints solid name and type."
    "\n  Note: all culling, if any, is switched off so all volumes are seen."
    "\n  Now printing with verbosity " << verbosity << G4endl;
}

void G4ASCIITreeSceneHandler::EndModeling () {
  fLVSet.clear();
  fReplicaSet.clear();
  G4cout << "G4ASCIITreeSceneHandler::EndModeling" << G4endl;
  G4VTreeSceneHandler::EndModeling ();  // To re-use "culling off" code.
}

void G4ASCIITreeSceneHandler::Dump (const G4VSolid& solid) {

  const G4ASCIITree* pSystem = (G4ASCIITree*)GetGraphicsSystem();
  const G4int verbosity = pSystem->GetVerbosity();
  const G4int detail = verbosity % 10;

  if (verbosity < 10 && fReplicaSet.find(fpCurrentPV) != fReplicaSet.end()) {
    // Ignore if an already treated replica.  (Assumes that the model
    // which has invoked this function is a G4PhysicalVolumeModel - we
    // check this by testing fpCurrentPV.)
    if (fpCurrentPV) {
      ((G4PhysicalVolumeModel*)fpModel)->CurtailDescent();
      return;
    }
  }

  // Print indented text...
  for (G4int i = 0; i < fCurrentDepth; i++ ) G4cout << "  ";
  G4cout << "\"" << fpCurrentPV->GetName()
	 << "\", copy no. " << fpCurrentPV->GetCopyNo();
  if (detail >= 1) {
    G4cout << ", belongs to logical volume \""
	   << fpCurrentLV->GetName() << "\"";
  }
  if (detail >= 2) {
    G4cout << " and is composed of solid \""
           << fpCurrentLV->GetSolid()->GetName()
	   << "\" of type \""
	   << fpCurrentLV->GetSolid()->GetEntityType() << "\"";
  }

  if (fpCurrentPV->IsReplicated()) {
    fReplicaSet.insert(fpCurrentPV);  // Record new replica volume.
    if (verbosity < 10) {
      // Add printing for replicas (when replicas are ignored)...
      EAxis axis;
      G4int nReplicas;
      G4double width;
      G4double offset;
      G4bool consuming;
      fpCurrentPV->GetReplicationData(axis,nReplicas,width,offset,consuming);
      G4VPVParameterisation* pP = fpCurrentPV->GetParameterisation();
      G4cout << " (" << nReplicas;
      if (pP) {
	G4cout << " parametrised volumes)";
      }
      else {
	G4cout << " replicas)";
      }
    }
  }
  else {
    if (fLVSet.find(fpCurrentLV) != fLVSet.end()) {
      if (verbosity <  10) {
	// Add printing for repeated logical volume...
	G4cout << " (repeated logical volume)";
	// Ignore if an already treated logical volume.
	if (fpCurrentPV) {
	  ((G4PhysicalVolumeModel*)fpModel)->CurtailDescent();
	  G4cout << G4endl;
	  return;
	}
      }
    }
  }

  if (fLVSet.find(fpCurrentLV) == fLVSet.end()) {
    fLVSet.insert(fpCurrentLV);  // Record new logical volume.
  }

  G4cout << G4endl;
  return;
}

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
// $Id: G4ASCIITreeSceneHandler.cc,v 1.32 2006/11/05 20:44:28 allison Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// 
// John Allison  5th April 2001
// A scene handler to dump geometry hierarchy.
// Based on a provisional G4ASCIITreeGraphicsScene (was in modeling).

#include "G4ASCIITreeSceneHandler.hh"

#include "G4ASCIITree.hh"
#include "G4ASCIITreeMessenger.hh"
#include "G4VSolid.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4Polyhedron.hh"
#include "G4UnitsTable.hh"
#include "G4Material.hh"
#include "G4Scene.hh"
#include "G4ModelingParameters.hh"
#include "G4PhysicalVolumeMassScene.hh"
#include "G4VSensitiveDetector.hh"
#include "G4VReadOutGeometry.hh"
#include "G4TransportationManager.hh"

G4ASCIITreeSceneHandler::G4ASCIITreeSceneHandler
(G4VGraphicsSystem& system,
 const G4String& name):
  G4VTreeSceneHandler(system, name),
  fpLastPV (0),
  fPVPCount (0),
  fpOutFile (0)
{}

G4ASCIITreeSceneHandler::~G4ASCIITreeSceneHandler () {}

void G4ASCIITreeSceneHandler::BeginModeling () {

  G4VTreeSceneHandler::BeginModeling ();  // To re-use "culling off" code.

  const G4ASCIITree* pSystem = (G4ASCIITree*)GetGraphicsSystem();
  const G4String& outFileName = pSystem -> GetOutFileName();
  if (outFileName == "G4cout") {
    fpOutFile = &G4cout;
  } else {
    fOutFile.open (outFileName);
    fpOutFile = &fOutFile;
  }

  G4cout << "G4ASCIITreeSceneHandler::BeginModeling: writing to ";
  if (outFileName == "G4cout") {
    G4cout << "G4 standard output (G4cout)";
  } else {
    G4cout << "file \"" << outFileName << "\"";
  }
  G4cout << G4endl;

  WriteHeader (G4cout); G4cout << G4endl;
  if (outFileName != "G4cout") {
    WriteHeader (fOutFile); fOutFile << std::endl;
  }
}

void G4ASCIITreeSceneHandler::WriteHeader (std::ostream& os)
{
  const G4ASCIITree* pSystem = (G4ASCIITree*)GetGraphicsSystem();
  const G4int verbosity = pSystem->GetVerbosity();
  const G4int detail = verbosity % 10;
  os << "#  Set verbosity with \"/vis/ASCIITree/verbose <verbosity>\":";
  for (size_t i = 0;
       i < G4ASCIITreeMessenger::fVerbosityGuidance.size(); ++i) {
    os << "\n#  " << G4ASCIITreeMessenger::fVerbosityGuidance[i];
  }
  os << "\n#  Now printing with verbosity " << verbosity;
  os << "\n#  Format is: PV:n";
  if (detail >= 1) os << " / LV (SD,RO)";
  if (detail >= 2) os << " / Solid(type)";
  if (detail >= 3) os << ", volume, density";
  if (detail >= 5) os << ", daughter-subtracted volume and mass";
  os <<
    "\n#  Abbreviations: PV = Physical Volume,     LV = Logical Volume,"
    "\n#                 SD = Sensitive Detector,  RO = Read Out Geometry.";
}

void G4ASCIITreeSceneHandler::EndModeling () {
  const G4ASCIITree* pSystem = (G4ASCIITree*) GetGraphicsSystem();
  const G4int verbosity = pSystem->GetVerbosity();
  const G4int detail = verbosity % 10;
  const G4String& outFileName = pSystem -> GetOutFileName();

  if (detail >= 4) {
    G4cout << "Calculating mass(es)..." << G4endl;
    const std::vector<G4VModel*>& models = fpScene->GetRunDurationModelList();
    std::vector<G4VModel*>::const_iterator i;
    for (i = models.begin(); i != models.end(); ++i) {
      G4PhysicalVolumeModel* pvModel =
	dynamic_cast<G4PhysicalVolumeModel*>(*i);
      if (pvModel) {
	if (pvModel->GetTopPhysicalVolume() ==
	    G4TransportationManager::GetTransportationManager()
	    ->GetNavigatorForTracking()->GetWorldVolume()) {
	  const G4ModelingParameters* tempMP =
	    pvModel->GetModelingParameters();
	  G4ModelingParameters mp;  // Default - no culling.
	  pvModel->SetModelingParameters (&mp);
	  G4PhysicalVolumeMassScene massScene(pvModel);
	  pvModel->DescribeYourselfTo (massScene);
	  G4double volume = massScene.GetVolume();
	  G4double mass = massScene.GetMass();

	  G4cout << "Overall volume of \""
		 << pvModel->GetTopPhysicalVolume()->GetName()
		 << "\":"
		 << pvModel->GetTopPhysicalVolume()->GetCopyNo()
		 << ", is "
		 << G4BestUnit (volume, "Volume")
		 << " and the daughter-included mass";
	  G4int requestedDepth = pvModel->GetRequestedDepth();
	  if (requestedDepth == G4PhysicalVolumeModel::UNLIMITED) {
	    G4cout << " to unlimited depth";
	  } else {
	    G4cout << ", ignoring daughters at depth "
		   << requestedDepth
		   << " and below,";
	  }
	  G4cout << " is " << G4BestUnit (mass, "Mass")
		 << G4endl;

	  pvModel->SetModelingParameters (tempMP);
	}
      }
    }
  }

  if (outFileName != "G4cout") {
    fOutFile.close();
    G4cout << "Output file \"" << outFileName << "\" closed." << G4endl;
  }
  fpLastPV = 0;
  fPVPCount = 0;
  fLVSet.clear();
  fReplicaSet.clear();
  G4cout << "G4ASCIITreeSceneHandler::EndModeling" << G4endl;
  G4VTreeSceneHandler::EndModeling ();  // To re-use "culling off" code.
}

void G4ASCIITreeSceneHandler::RequestPrimitives(const G4VSolid& solid) {
  
  G4PhysicalVolumeModel* pPVModel =
    dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
  if (!pPVModel) return;  // Not from a G4PhysicalVolumeModel.

  // This call comes from a G4PhysicalVolumeModel.  drawnPVPath is
  // the path of the current drawn (non-culled) volume in terms of
  // drawn (non-culled) ancesters.  Each node is identified by a
  // PVNodeID object, which is a physical volume and copy number.  It
  // is a vector of PVNodeIDs corresponding to the geometry hierarchy
  // actually selected, i.e., not culled.
  typedef G4PhysicalVolumeModel::G4PhysicalVolumeNodeID PVNodeID;
  typedef std::vector<PVNodeID> PVPath;
  const PVPath& drawnPVPath = pPVModel->GetDrawnPVPath();
  //G4int currentDepth = pPVModel->GetCurrentDepth();
  G4VPhysicalVolume* pCurrentPV = pPVModel->GetCurrentPV();
  G4LogicalVolume* pCurrentLV = pPVModel->GetCurrentLV();
  G4Material* pCurrentMaterial = pPVModel->GetCurrentMaterial();

  G4ASCIITree* pSystem = (G4ASCIITree*)GetGraphicsSystem();
  G4int verbosity = pSystem->GetVerbosity();
  G4int detail = verbosity % 10;
  const G4String& outFileName = pSystem -> GetOutFileName();

  if (pCurrentPV != fpLastPV) {
    fpLastPV = pCurrentPV;
    fPVPCount = 0;
  }

  if (verbosity < 10 && pCurrentPV->IsReplicated()) {
    // See if this has been a replica found before with same mother LV...
    PVPath::const_reverse_iterator thisID = drawnPVPath.rbegin();
    PVPath::const_reverse_iterator motherID = ++drawnPVPath.rbegin();
    G4bool ignore = false;
    for (ReplicaSetIterator i = fReplicaSet.begin(); i != fReplicaSet.end();
	 ++i) {
      if (i->back().GetPhysicalVolume()->GetLogicalVolume() ==
	  thisID->GetPhysicalVolume()->GetLogicalVolume()) {
	// For each one previously found (if more than one, they must
	// have different mothers)...
	// To avoid compilation errors on VC++ .Net 7.1...
	// Previously:
	//   PVNodeID previousMotherID = ++(i->rbegin());
	// (Should that have been: PVNodeID::const_iterator previousMotherID?)
	// Replace
	//   previousMotherID == i->rend()
	// by
	//   i->size() <= 1
	// Replace
	//   previousMotherID != i->rend()
	// by
	//   i->size() > 1
	// Replace
	//   previousMotherID->
	// by
	//   (*i)[i->size() - 2].
	if (motherID == drawnPVPath.rend() &&
	    i->size() <= 1)
	  ignore = true;  // Both have no mother.
	if (motherID != drawnPVPath.rend() &&
	    i->size() > 1 &&
	    motherID->GetPhysicalVolume()->GetLogicalVolume() ==
	    (*i)[i->size() - 2].GetPhysicalVolume()->GetLogicalVolume())
	  ignore = true;  // Same mother LV...
      }
    }
    if (ignore) {
      pPVModel->CurtailDescent();
      return;
    }
  }

  // Print indented text...
  if (pCurrentPV) {
    for (size_t i = 0; i < drawnPVPath.size(); i++ ) *fpOutFile << "  ";

    *fpOutFile << "\"" << pCurrentPV->GetName()
	       << "\":" << pCurrentPV->GetCopyNo();

    if (pCurrentPV->IsReplicated()) {
      if (verbosity < 10) {
	// Add printing for replicas (when replicas are ignored)...
	EAxis axis;
	G4int nReplicas;
	G4double width;
	G4double offset;
	G4bool consuming;
	pCurrentPV->GetReplicationData(axis,nReplicas,width,offset,consuming);
	G4VPVParameterisation* pP = pCurrentPV->GetParameterisation();
	if (pP) {
	  if (detail < 3) {
	    fReplicaSet.insert(drawnPVPath);
	    *fpOutFile << " (" << nReplicas << " parametrised volumes)";
	  }
	}
	else {
	  fReplicaSet.insert(drawnPVPath);
	  *fpOutFile << " (" << nReplicas << " replicas)";
	}
      }
    }
    else {
      if (fLVSet.find(pCurrentLV) != fLVSet.end()) {
	if (verbosity <  10) {
	  // Add printing for repeated LV...
	  *fpOutFile << " (repeated LV)";
	  // ...and curtail descent.
	  pPVModel->CurtailDescent();
	}
      }
    }

    if (detail >= 1) {
      *fpOutFile << " / \""
		 << pCurrentLV->GetName() << "\"";
      G4VSensitiveDetector* sd = pCurrentLV->GetSensitiveDetector();
      if (sd) {
	*fpOutFile << " (SD=\"" << sd->GetName() << "\"";
	G4VReadOutGeometry* roGeom = sd->GetROgeometry();
	if (roGeom) {
	  *fpOutFile << ",RO=\"" << roGeom->GetName() << "\"";
	}
	*fpOutFile << ")";
      }
    }
  }  // if (pCurrentPV)

  if (detail >= 2) {
    *fpOutFile << " / \""
	       << solid.GetName()
	       << "\"("
	       << solid.GetEntityType() << ")";
  }

  if (pCurrentPV) {
    if (detail >= 3) {
      *fpOutFile << ", "
		 << G4BestUnit(((G4VSolid&)solid).GetCubicVolume(),"Volume")
		 << ", ";
      if (pCurrentMaterial) {
	*fpOutFile
	  << G4BestUnit(pCurrentMaterial->GetDensity(), "Volumic Mass")
		   << " (" << pCurrentMaterial->GetName() << ")";
      } else {
	*fpOutFile << "(No material)";
      }
    }

    if (detail >= 5) {
      if (pCurrentMaterial) {
	G4Material* pMaterial = const_cast<G4Material*>(pCurrentMaterial);
	// ...and find daughter-subtracted mass...
	G4double daughter_subtracted_mass = pCurrentLV->GetMass
	  (pCurrentPV->IsParameterised(),  // Force if parametrised.
	   false,  // Do not propagate - calculate for this volume minus
	           // volume of daughters.
	   pMaterial);
	G4double daughter_subtracted_volume =
	  daughter_subtracted_mass / pCurrentMaterial->GetDensity();
	*fpOutFile << ", "
		   << G4BestUnit(daughter_subtracted_volume,"Volume")
		   << ", "
		   << G4BestUnit(daughter_subtracted_mass,"Mass");
      }
    }

    if (fLVSet.find(pCurrentLV) == fLVSet.end()) {
      fLVSet.insert(pCurrentLV);  // Record new logical volume.
    }
  }  // if (pCurrentPV)

  if (outFileName == "G4cout") {
    G4cout << G4endl;
  } else {
    *fpOutFile << std::endl;
  }

  return;
}

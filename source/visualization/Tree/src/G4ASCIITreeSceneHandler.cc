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
// $Id: G4ASCIITreeSceneHandler.cc,v 1.16 2004/11/11 16:03:15 johna Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
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
#include "G4Polyhedron.hh"
#include "G4UnitsTable.hh"
#include "G4Material.hh"
#include "G4Scene.hh"
#include "G4ModelingParameters.hh"
#include "G4PhysicalVolumeMassScene.hh"

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
  os <<
    "#  Set verbosity with \"/vis/ASCIITree/verbose <verbosity>\":"
    "\n#  <  10: - does not print daughters of repeated placements."
    "\n#         - does not repeat replicas."
    "\n#  >= 10: prints all physical volumes."
    "\n#  The level of detail is given by the units (verbosity%10):"
    "\n#  >=  0: prints physical volume name."
    "\n#  >=  1: prints logical volume name."
    "\n#  >=  2: prints solid name and type."
    "\n#  >=  3: prints volume and density of solid."
    "\n#  >=  4: calculates and prints mass(es) of volume(s) in scene."
    "\n#  Note: by default, culling is switched off so all volumes are seen.";
  if (detail >=4) {
    os <<
      "\n#  Note: the mass is calculated for each physical volume in the scene,"
      "\n#  taking into account daughters up to the depth specified.  Culling "
      "\n#  is ignored.  If you want the mass of a particular subtree:"
      "\n#    /vis/drawVolume <name-of-physical-volume-at-top-of-subtree>"
      "\n#    /vis/viewer/flush";
  }

  os << "\n#  Now printing with verbosity " << verbosity;
  os << "\n#  Format is: PV:n";
  if (detail >= 1) os << " / LV";
  if (detail >= 2) os << " / Solid(type)";
  if (detail >= 3) os << ", volume, density";
  os << "\n#  where PV = Physical Volume";
  if (detail <2) os << " and"; else os << ",";
  os << " n = copy number";
  if (detail >= 2) os << " and LV = Logical Volume";
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
    G4PhysicalVolumeModel* pvModel;
    for (i = models.begin(); i != models.end(); ++i) {
      if ((pvModel = (*i)->GetG4PhysicalVolumeModel())) {
	const G4ModelingParameters* tempMP = pvModel->GetModelingParameters();
	G4ModelingParameters mp;  // Default - no culling.
	pvModel->SetModelingParameters (&mp);
	G4PhysicalVolumeMassScene massScene;
	massScene.EstablishSpecials (*pvModel);
	pvModel->DescribeYourselfTo (massScene);
	G4double volume = massScene.GetVolume();
	G4double mass = massScene.GetMass();

	G4cout << "Overall volume of \""
	       << pvModel->GetTopPhysicalVolume()->GetName()
	       << "\":"
	       << pvModel->GetTopPhysicalVolume()->GetCopyNo()
	       << ", is "
	       << G4BestUnit (volume, "Volume")
	       << "\nMass of tree to ";
	G4int requestedDepth = pvModel->GetRequestedDepth();
	if (requestedDepth == G4PhysicalVolumeModel::UNLIMITED) {
	  G4cout << "unlimited depth";
	} else {
	  G4cout << "depth " << requestedDepth;
	}
	G4cout << " is " << G4BestUnit (mass, "Mass")
	       << G4endl;

	pvModel->SetModelingParameters (tempMP);
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
  
  const G4ASCIITree* pSystem = (G4ASCIITree*)GetGraphicsSystem();
  const G4int verbosity = pSystem->GetVerbosity();
  const G4int detail = verbosity % 10;
  const G4String& outFileName = pSystem -> GetOutFileName();

  if (fpCurrentPV != fpLastPV) {
    fpLastPV = fpCurrentPV;
    fPVPCount = 0;
  }

  if (verbosity < 10 && fReplicaSet.find(fpCurrentPV) != fReplicaSet.end()) {
    // Ignore if an already treated replica.
    G4PhysicalVolumeModel* pPVM = fpModel->GetG4PhysicalVolumeModel();
    if (pPVM) {
      pPVM->CurtailDescent();
      return;
    } else {
      G4Exception
	("G4ASCIITreeSceneHandler::RequestPrimitives:"
	 "not a physical volume model");  // Shouldn't happen.
    }
  }

  // Print indented text...
  for (G4int i = 0; i < fCurrentDepth; i++ ) *fpOutFile << "  ";

  *fpOutFile << "\"" << fpCurrentPV->GetName()
	     << "\":" << fpCurrentPV->GetCopyNo();

  if (fpCurrentPV->IsReplicated()) {
    if (verbosity < 10) {
      // Add printing for replicas (when replicas are ignored)...
      EAxis axis;
      G4int nReplicas;
      G4double width;
      G4double offset;
      G4bool consuming;
      fpCurrentPV->GetReplicationData(axis,nReplicas,width,offset,consuming);
      G4VPVParameterisation* pP = fpCurrentPV->GetParameterisation();
      if (pP) {
	if (detail < 3) {
	  fReplicaSet.insert(fpCurrentPV);
	  *fpOutFile << " (" << nReplicas << " parametrised volumes)";
	}
      }
      else {
	fReplicaSet.insert(fpCurrentPV);
	*fpOutFile << " (" << nReplicas << " replicas)";
      }
    }
  }
  else {
    if (fLVSet.find(fpCurrentLV) != fLVSet.end()) {
      if (verbosity <  10) {
	// Add printing for repeated placement...
	*fpOutFile << " (repeated placement)";
	// ...and curtail descent.
	((G4PhysicalVolumeModel*)fpModel)->CurtailDescent();
      }
    }
  }

  if (detail >= 1) {
    *fpOutFile << " / \""
	       << fpCurrentLV->GetName() << "\"";
  }

  if (detail >= 2) {
    *fpOutFile << " / \""
	       << solid.GetName()
	       << "\"("
	       << solid.GetEntityType() << ")";
  }

  if (detail >= 3) {
    G4Polyhedron* pPolyhedron = solid.GetPolyhedron();
    if (pPolyhedron) {
      G4Material* pMaterial;
      G4VPVParameterisation* pP = fpCurrentPV->GetParameterisation();
      if (pP) {
        pMaterial = pP -> ComputeMaterial (fPVPCount++, fpCurrentPV);
      } else {
	pMaterial = fpCurrentLV->GetMaterial();
      }
      *fpOutFile << ", "
		 << G4BestUnit(pPolyhedron->GetVolume(),"Volume")
		 << ", "
		 << G4BestUnit(pMaterial->GetDensity(), "Volumic Mass");
    } else {
      *fpOutFile << " (volume not available)";
    }
  }

  if (fLVSet.find(fpCurrentLV) == fLVSet.end()) {
    fLVSet.insert(fpCurrentLV);  // Record new logical volume.
  }

  if (outFileName == "G4cout") {
    G4cout << G4endl;
  } else {
    *fpOutFile << std::endl;
  }

  return;
}

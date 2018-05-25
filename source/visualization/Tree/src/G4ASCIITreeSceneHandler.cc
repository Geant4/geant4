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
// $Id: G4ASCIITreeSceneHandler.cc 110130 2018-05-16 06:43:10Z gcosmo $
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
#include "G4AttCheck.hh"
#include "G4AttValue.hh"

G4ASCIITreeSceneHandler::G4ASCIITreeSceneHandler
(G4VGraphicsSystem& system,
 const G4String& name):
  G4VTreeSceneHandler(system, name),
  fpOutFile(0),
  fpLastPV(0),
  fLastCopyNo(-99),
  fLastNonSequentialCopyNo(-99)
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

  static G4bool firstTime = true;
  if (firstTime) {
    firstTime = false;
    G4cout << "G4ASCIITreeSceneHandler::BeginModeling: writing to ";
    if (outFileName == "G4cout") {
      G4cout << "G4 standard output (G4cout)";
    } else {
      G4cout << "file \"" << outFileName << "\"";
    }
    G4cout << G4endl;

    WriteHeader (G4cout); G4cout << G4endl;
  }
  
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
  if (detail >= 6) os << ", physical volume dump";
  if (detail >= 7) os << ", polyhedron dump";
  os <<
    "\n#  Abbreviations: PV = Physical Volume,     LV = Logical Volume,"
    "\n#                 SD = Sensitive Detector,  RO = Read Out Geometry.";
}

void G4ASCIITreeSceneHandler::EndModeling () {
  const G4ASCIITree* pSystem = (G4ASCIITree*) GetGraphicsSystem();
  const G4int verbosity = pSystem->GetVerbosity();
  const G4int detail = verbosity % 10;
  const G4String& outFileName = pSystem -> GetOutFileName();

  // Output left over copy number, if any...
  if (fLastCopyNo != fLastNonSequentialCopyNo) {
    if (fLastCopyNo == fLastNonSequentialCopyNo + 1) *fpOutFile << ',';
    else *fpOutFile << '-';
    *fpOutFile << fLastCopyNo;
  }
  // Output outstanding rest of line, if any...
  if (!fRestOfLine.str().empty()) *fpOutFile << fRestOfLine.str();
  fRestOfLine.str("");
  fpLastPV = 0;
  fLastPVName.clear();
  fLastCopyNo = -99;
  fLastNonSequentialCopyNo = -99;

  // This detail to G4cout regardless of outFileName...
  if (detail >= 4) {
    G4cout << "Calculating mass(es)..." << G4endl;
    const std::vector<G4Scene::Model>& models = fpScene->GetRunDurationModelList();
    std::vector<G4Scene::Model>::const_iterator i;
    for (i = models.begin(); i != models.end(); ++i) {
      G4PhysicalVolumeModel* pvModel =
	dynamic_cast<G4PhysicalVolumeModel*>(i->fpModel);
      if (pvModel) {
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

  if (outFileName != "G4cout") {
    fOutFile.close();
    G4cout << "Output file \"" << outFileName << "\" closed." << G4endl;
  }
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
  // The following typedef's already set in header file...
  //typedef G4PhysicalVolumeModel::G4PhysicalVolumeNodeID PVNodeID;
  //typedef std::vector<PVNodeID> PVPath;
  const PVPath& drawnPVPath = pPVModel->GetDrawnPVPath();
  //G4int currentDepth = pPVModel->GetCurrentDepth();
  G4VPhysicalVolume* pCurrentPV = pPVModel->GetCurrentPV();
  const G4String& currentPVName = pCurrentPV->GetName();
  const G4int currentCopyNo = pCurrentPV->GetCopyNo();
  G4LogicalVolume* pCurrentLV = pPVModel->GetCurrentLV();
  G4Material* pCurrentMaterial = pPVModel->GetCurrentMaterial();

  G4ASCIITree* pSystem = (G4ASCIITree*)GetGraphicsSystem();
  G4int verbosity = pSystem->GetVerbosity();
  G4int detail = verbosity % 10;

  // If verbosity < 10 suppress unnecessary repeated printing.
  // Repeated simple replicas can always be suppressed.
  // Paramaterisations can only be suppressed if verbosity < 3, since their
  // size, density, etc., are in principle different.
  const G4bool isParameterised = pCurrentPV->GetParameterisation();
  const G4bool isSimpleReplica = pCurrentPV->IsReplicated() && !isParameterised;
  const G4bool isAmenableToSupression =
  (verbosity < 10 && isSimpleReplica) || (verbosity < 3 && isParameterised);
  if (isAmenableToSupression) {
    // See if this has been found before with same mother LV...
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

  // Now suppress printing for volumes with the same name but different
  // copy number - but not if they are parameterisations that have not been
  // taken out above (those that are "amenable to suppression" will have been
  // taken out).
  if (verbosity < 10 && !isParameterised &&
      currentPVName == fLastPVName &&
      currentCopyNo != fLastCopyNo) {
    // Check...
    if (isAmenableToSupression) {
      G4Exception("G4ASCIITreeSceneHandler::RequestPrimitives",
		  "vistree0001",
		  JustWarning,
		  "Volume amenable to suppressed printing unexpected");
    }
    // Check mothers are identical...
    else if (pCurrentLV == (fpLastPV? fpLastPV->GetLogicalVolume(): 0)) {
      if (currentCopyNo != fLastCopyNo + 1) {
	// Non-sequential copy number...
	*fpOutFile << ',' << currentCopyNo;
	fLastNonSequentialCopyNo = currentCopyNo;
      }
      fLastCopyNo = currentCopyNo;
      pPVModel->CurtailDescent();
      return;
    }
  }
  fpLastPV = pCurrentPV;

  // High verbosity or a new or otherwise non-amenable volume...
  // Output copy number, if any, from previous invocation...
  if (fLastCopyNo != fLastNonSequentialCopyNo) {
    if (fLastCopyNo == fLastNonSequentialCopyNo + 1) *fpOutFile << ',';
    else *fpOutFile << '-';
    *fpOutFile << fLastCopyNo;
  }
  // Output rest of line, if any, from previous invocation...
  if (!fRestOfLine.str().empty()) *fpOutFile << fRestOfLine.str();
  fRestOfLine.str("");
  fLastPVName = currentPVName;
  fLastCopyNo = currentCopyNo;
  fLastNonSequentialCopyNo = currentCopyNo;
  // Start next line...
  // Indent according to drawn path depth...
  for (size_t i = 0; i < drawnPVPath.size(); i++ ) *fpOutFile << "  ";
  *fpOutFile << "\"" << currentPVName
	     << "\":" << currentCopyNo;

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
	  if (nReplicas > 2) fRestOfLine << '-';
	  else fRestOfLine << ',';
	  fRestOfLine << nReplicas - 1
		      << " (" << nReplicas << " parametrised volumes)";
	}
      }
      else {
	fReplicaSet.insert(drawnPVPath);
	if (nReplicas > 2) fRestOfLine << '-';
	else fRestOfLine << ',';
	fRestOfLine << nReplicas - 1
		    << " (" << nReplicas << " replicas)";
      }
    }
  } else {
    if (fLVSet.find(pCurrentLV) != fLVSet.end()) {
      if (verbosity <  10) {
	// Add printing for repeated LV (if it has daughters)...
	if (pCurrentLV->GetNoDaughters()) fRestOfLine << " (repeated LV)";
	// ...and curtail descent.
	pPVModel->CurtailDescent();
      }
    }
  }

  if (detail >= 1) {
    fRestOfLine << " / \""
		<< pCurrentLV->GetName() << "\"";
    G4VSensitiveDetector* sd = pCurrentLV->GetSensitiveDetector();
    if (sd) {
      fRestOfLine << " (SD=\"" << sd->GetName() << "\"";
      G4VReadOutGeometry* roGeom = sd->GetROgeometry();
      if (roGeom) {
	fRestOfLine << ",RO=\"" << roGeom->GetName() << "\"";
      }
      fRestOfLine << ")";
    }
  }

  if (detail >= 2) {
    fRestOfLine << " / \""
	       << solid.GetName()
	       << "\"("
	       << solid.GetEntityType() << ")";
  }

  if (detail >= 3) {
    fRestOfLine << ", "
		<< G4BestUnit(((G4VSolid&)solid).GetCubicVolume(),"Volume")
		<< ", ";
    if (pCurrentMaterial) {
      fRestOfLine
	<< G4BestUnit(pCurrentMaterial->GetDensity(), "Volumic Mass")
	<< " (" << pCurrentMaterial->GetName() << ")";
    } else {
      fRestOfLine << "(No material)";
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
      fRestOfLine << ", "
		  << G4BestUnit(daughter_subtracted_volume,"Volume")
		  << ", "
		  << G4BestUnit(daughter_subtracted_mass,"Mass");
    }
  }

  if (detail >= 6) {
    std::vector<G4AttValue>* attValues = pPVModel->CreateCurrentAttValues();
    const std::map<G4String,G4AttDef>* attDefs = pPVModel->GetAttDefs();
    fRestOfLine << '\n' << G4AttCheck(attValues, attDefs);
    delete attValues;
  }

  if (detail >= 7) {
    G4Polyhedron* polyhedron = solid.GetPolyhedron();
    fRestOfLine << "\nLocal polyhedron coordinates:\n" << *polyhedron;
    G4Transform3D* transform = pPVModel->GetCurrentTransform();
    polyhedron->Transform(*transform);
    fRestOfLine << "\nGlobal polyhedron coordinates:\n" << *polyhedron;
  }

  if (fLVSet.find(pCurrentLV) == fLVSet.end()) {
    fLVSet.insert(pCurrentLV);  // Record new logical volume.
  }

  fRestOfLine << std::endl;

  return;
}

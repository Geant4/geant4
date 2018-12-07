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

// /vis/touchable commands - John Allison  14th May 2014

#include "G4VisCommandsTouchable.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4TransportationManager.hh"
#include "G4TouchableUtils.hh"
#include "G4PhysicalVolumesSearchScene.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4AttCheck.hh"

G4VisCommandsTouchable::G4VisCommandsTouchable()
{
  fpCommandDump = new G4UIcmdWithoutParameter("/vis/touchable/dump",this);
  fpCommandDump->SetGuidance("Dump touchable attributes.");
  fpCommandDump->SetGuidance
  ("Use \"/vis/set/touchable\" to set current touchable.");
  fpCommandDump->SetGuidance
  ("Use \"/vis/touchable/set\" to set attributes.");

  fpCommandFindPath = new G4UIcmdWithAString("/vis/touchable/findPath",this);
  fpCommandFindPath->SetGuidance
  ("Prints the path to touchable and its logical volume mother"
   "\ngiven a physical volume name.");
  fpCommandFindPath -> SetGuidance
  ("A search of all worlds is made and all physical volume names are"
   "\nmatched against the argument of this command.  If this is of the"
   "\nform \"/regexp/\", where regexp is a regular expression (see C++ regex),"
   "\nthe physical volume name is matched against regexp by the usual rules"
   "\nof regular expression matching. Otherwise an exact match is required."
   "\nFor example, \"/Shap/\" matches \"Shape1\" and \"Shape2\".");
  fpCommandFindPath -> SetGuidance
  ("It may help to see a textual representation of the geometry hierarchy of"
   "\nthe worlds. Try \"/vis/drawTree [worlds]\" or one of the driver/browser"
   "\ncombinations that have the required functionality, e.g., HepRep.");
}

G4VisCommandsTouchable::~G4VisCommandsTouchable() {
  delete fpCommandDump;
}

G4String G4VisCommandsTouchable::GetCurrentValue(G4UIcommand*) {
  return "";
}

void G4VisCommandsTouchable::SetNewValue
(G4UIcommand* command, G4String newValue)
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

  if (command == fpCommandDump) {

    G4PhysicalVolumeModel::TouchableProperties properties =
    G4TouchableUtils::FindTouchableProperties(fCurrentTouchableProperties.fTouchablePath);
    if (properties.fpTouchablePV) {
      G4PhysicalVolumeModel pvModel
      (properties.fpTouchablePV,
       G4PhysicalVolumeModel::UNLIMITED,
       properties.fTouchableGlobalTransform,
       nullptr, // Modelling parameters (not used)
       true, // use full extent (prevents calculating own extent, which crashes)
       properties.fTouchableBaseFullPVPath);
      const std::map<G4String,G4AttDef>* attDefs = pvModel.GetAttDefs();
      std::vector<G4AttValue>* attValues = pvModel.CreateCurrentAttValues();
      G4cout << G4AttCheck(attValues,attDefs);
      delete attValues;
      G4Polyhedron* polyhedron =
      properties.fpTouchablePV->GetLogicalVolume()->GetSolid()->GetPolyhedron();
      G4cout << "\nLocal polyhedron coordinates:\n" << *polyhedron;
      G4Transform3D* transform = pvModel.GetCurrentTransform();
      polyhedron->Transform(*transform);
      G4cout << "\nGlobal polyhedron coordinates:\n" << *polyhedron;
    } else {
      G4cout << "Touchable not found." << G4endl;
    }
    return;

  } else if (command == fpCommandFindPath) {

    std::vector<G4PhysicalVolumesSearchScene::Findings> findingsVector;
    std::vector<G4VPhysicalVolume*>::iterator iterWorld =
    transportationManager->GetWorldsIterator();
    for (size_t i = 0; i < nWorlds; ++i, ++iterWorld) {
      G4PhysicalVolumeModel searchModel (*iterWorld);  // Unlimited depth.
      G4ModelingParameters mp;  // Default - no culling.
      searchModel.SetModelingParameters (&mp);
      G4PhysicalVolumesSearchScene searchScene (&searchModel, newValue);
      searchModel.DescribeYourselfTo (searchScene);  // Initiate search.
      for (const auto& findings: searchScene.GetFindings()) {
        findingsVector.push_back(findings);
      }
    }
    for (const auto& findings: findingsVector) {
      G4cout
      <<  findings.fFoundBasePVPath
      << ' ' << findings.fpFoundPV->GetName()
      << ' ' <<  findings.fFoundPVCopyNo
      << " (mother logical volume: "
      << findings.fpFoundPV->GetMotherLogical()->GetName()
      << ')'
      << G4endl;
    }
    if (findingsVector.size()) {
      G4cout
      << "Use this to set a particular touchable with \"/vis/set/touchable <path>\""
      << "\nor to see overlaps: \"/vis/drawLogicalVolume <mother-logical-volume-name>\""
      << G4endl;
    } else {
      G4cout << newValue << " not found" << G4endl;
    }

  } else {

    if (verbosity >= G4VisManager::errors) {
      G4cerr <<
      "ERROR: G4VisCommandsTouchable::SetNewValue: unrecognised command."
      << G4endl;
    }
    return;
  }
}

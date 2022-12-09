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

#include "G4UImanager.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4TransportationManager.hh"
#include "G4TouchableUtils.hh"
#include "G4PhysicalVolumesSearchScene.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4AttCheck.hh"
#include "G4AxesModel.hh"

#define G4warn G4cout

G4VisCommandsTouchable::G4VisCommandsTouchable()
{
  G4bool omitable;

  fpCommandCentreAndZoomInOn = new G4UIcmdWithoutParameter("/vis/touchable/centreAndZoomInOn",this);
  fpCommandCentreAndZoomInOn->SetGuidance ("Centre and zoom in on the current touchable.");
  fpCommandCentreAndZoomInOn->SetGuidance
  ("Use \"/vis/set/touchable\" to set current touchable.");
  fpCommandCentreAndZoomInOn->SetGuidance
  ("You may also need \"/vis/touchable/findPath\".");
  fpCommandCentreAndZoomInOn->SetGuidance
  ("Use \"/vis/touchable/set\" to set attributes.");

  fpCommandCentreOn = new G4UIcmdWithoutParameter("/vis/touchable/centreOn",this);
  fpCommandCentreOn->SetGuidance ("Centre the view on the current touchable.");
  // Pick up additional guidance from /vis/viewer/centreAndZoomInOn
  CopyGuidanceFrom(fpCommandCentreAndZoomInOn,fpCommandCentreOn,1);

  fpCommandDraw = new G4UIcmdWithABool("/vis/touchable/draw",this);
  fpCommandDraw->SetGuidance("Draw touchable.");
  fpCommandDraw->SetGuidance
  ("If parameter == true, also draw extent as a white wireframe box.");
  // Pick up additional guidance from /vis/viewer/centreAndZoomInOn
  CopyGuidanceFrom(fpCommandCentreAndZoomInOn,fpCommandDraw,1);
  fpCommandDraw->SetParameterName("extent", omitable = true);
  fpCommandDraw->SetDefaultValue(false);

  fpCommandDump = new G4UIcmdWithoutParameter("/vis/touchable/dump",this);
  fpCommandDump->SetGuidance("Dump touchable attributes.");
  // Pick up additional guidance from /vis/viewer/centreAndZoomInOn
  CopyGuidanceFrom(fpCommandCentreAndZoomInOn,fpCommandDump,1);

  fpCommandExtentForField = new G4UIcmdWithABool("/vis/touchable/extentForField",this);
  fpCommandExtentForField->SetGuidance("Set extent for field.");
  fpCommandExtentForField->SetGuidance("If parameter == true, also draw.");
  // Pick up additional guidance from /vis/viewer/centreAndZoomInOn
  CopyGuidanceFrom(fpCommandCentreAndZoomInOn,fpCommandExtentForField,1);
  fpCommandExtentForField->SetParameterName("draw", omitable = true);
  fpCommandExtentForField->SetDefaultValue(false);

  fpCommandFindPath = new G4UIcommand("/vis/touchable/findPath",this);
  fpCommandFindPath->SetGuidance
  ("Prints the path to touchable and its logical volume mother"
   "\ngiven a physical volume name and copy no.");
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
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("physical-volume-name", 's', omitable = true);
  parameter -> SetDefaultValue ("world");
  fpCommandFindPath -> SetParameter (parameter);
  parameter = new G4UIparameter ("copy-no", 'i', omitable = true);
  parameter -> SetGuidance ("If negative, matches any copy no.");
  parameter -> SetDefaultValue (-1);
  fpCommandFindPath -> SetParameter (parameter);

  fpCommandLocalAxes = new G4UIcmdWithoutParameter("/vis/touchable/localAxes",this);
  fpCommandLocalAxes->SetGuidance("Draw local axes.");
  // Pick up additional guidance from /vis/viewer/centreAndZoomInOn
  CopyGuidanceFrom(fpCommandCentreAndZoomInOn,fpCommandLocalAxes,1);

  fpCommandShowExtent = new G4UIcmdWithABool("/vis/touchable/showExtent",this);
  fpCommandShowExtent->SetGuidance("Print extent of touchable.");
  fpCommandShowExtent->SetGuidance("If parameter == true, also draw.");
  // Pick up additional guidance from /vis/viewer/centreAndZoomInOn
  CopyGuidanceFrom(fpCommandCentreAndZoomInOn,fpCommandShowExtent,1);
  fpCommandShowExtent->SetParameterName("draw", omitable = true);
  fpCommandShowExtent->SetDefaultValue(false);

  fpCommandVolumeForField = new G4UIcmdWithABool("/vis/touchable/volumeForField",this);
  fpCommandVolumeForField->SetGuidance("Set volume for field.");
  fpCommandVolumeForField->SetGuidance("If parameter == true, also draw.");
  // Pick up additional guidance from /vis/viewer/centreAndZoomInOn
  CopyGuidanceFrom(fpCommandCentreAndZoomInOn,fpCommandVolumeForField,1);
  fpCommandVolumeForField->SetParameterName("draw", omitable = true);
  fpCommandVolumeForField->SetDefaultValue(false);
}

G4VisCommandsTouchable::~G4VisCommandsTouchable() {
  delete fpCommandVolumeForField;
  delete fpCommandShowExtent;
  delete fpCommandLocalAxes;
  delete fpCommandFindPath;
  delete fpCommandExtentForField;
  delete fpCommandDump;
  delete fpCommandDraw;
  delete fpCommandCentreAndZoomInOn;
  delete fpCommandCentreOn;
}

G4String G4VisCommandsTouchable::GetCurrentValue(G4UIcommand*) {
  return "";
}

void G4VisCommandsTouchable::SetNewValue
(G4UIcommand* command, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn = verbosity >= G4VisManager::warnings;

  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  G4TransportationManager* transportationManager =
  G4TransportationManager::GetTransportationManager ();

  size_t nWorlds = transportationManager->GetNoWorlds();

  G4VPhysicalVolume* world = *(transportationManager->GetWorldsIterator());
  if (!world) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<
      "ERROR: G4VisCommandsTouchable::SetNewValue:"
      "\n  No world.  Maybe the geometry has not yet been defined."
      "\n  Try \"/run/initialize\""
      << G4endl;
    }
    return;
  }

  G4VViewer* currentViewer = fpVisManager -> GetCurrentViewer ();
  if (!currentViewer) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<
      "ERROR: No current viewer - \"/vis/viewer/list\" to see possibilities."
      << G4endl;
    }
    return;
  }

  G4Scene* currentScene = fpVisManager->GetCurrentScene();
  if (!currentScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<
      "ERROR: No current scene - \"/vis/scene/list\" to see possibilities."
      << G4endl;
    }
    return;
  }

  if (command == fpCommandCentreOn || command == fpCommandCentreAndZoomInOn) {

    // For twinkling...
    std::vector<std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>> touchables;

    G4PhysicalVolumeModel::TouchableProperties properties =
    G4TouchableUtils::FindTouchableProperties(fCurrentTouchableProperties.fTouchablePath);
    if (properties.fpTouchablePV) {
      // To handle parameterisations, set copy number
      properties.fpTouchablePV->SetCopyNo(properties.fCopyNo);
      G4PhysicalVolumeModel tempPVModel
      (properties.fpTouchablePV,
       G4PhysicalVolumeModel::UNLIMITED,
       properties.fTouchableGlobalTransform,
       nullptr, // Modelling parameters (not used)
       true, // use full extent (prevents calculating own extent, which crashes)
       properties.fTouchableBaseFullPVPath);
      touchables.push_back(properties.fTouchableFullPVPath);  // Only one in this case
      // Use a temporary scene in order to find vis extent
      G4Scene tempScene("Centre Scene");
      G4bool successful = tempScene.AddRunDurationModel(&tempPVModel,warn);
      if (!successful) return;
      if (verbosity >= G4VisManager::parameters) {
        G4cout
        << "Touchable " << fCurrentTouchableProperties.fTouchablePath
        << ",\n  has been added to temporary scene \"" << tempScene.GetName() << "\"."
        << G4endl;
      }

      const G4VisExtent& newExtent = tempScene.GetExtent();
      const G4ThreeVector& newTargetPoint = newExtent.GetExtentCentre();
      G4ViewParameters saveVP = currentViewer->GetViewParameters();
      G4ViewParameters newVP = saveVP;
      if (command == fpCommandCentreAndZoomInOn) {
        // Calculate the new zoom factor
        const G4double zoomFactor
        = currentScene->GetExtent().GetExtentRadius()/newExtent.GetExtentRadius();
        newVP.SetZoomFactor(zoomFactor);
      }
      // Change the target point
      const G4Point3D& standardTargetPoint = currentScene->GetStandardTargetPoint();
      newVP.SetCurrentTargetPoint(newTargetPoint - standardTargetPoint);

      // Interpolate
      auto keepVisVerbose = fpVisManager->GetVerbosity();
      fpVisManager->SetVerboseLevel(G4VisManager::errors);
      if (newVP != saveVP) InterpolateToNewView(currentViewer, saveVP, newVP);
      // ...and twinkle
      Twinkle(currentViewer,newVP,touchables);
      fpVisManager->SetVerboseLevel(keepVisVerbose);

      if (verbosity >= G4VisManager::confirmations) {
        G4cout
        << "Viewer \"" << currentViewer->GetName()
        << "\" centred ";
        if (fpCommandCentreAndZoomInOn) {
          G4cout << "and zoomed in";
        }
        G4cout << " on touchable\n" << fCurrentTouchableProperties.fTouchablePath
        << G4endl;
      }
      SetViewParameters(currentViewer, newVP);
    } else {
      G4warn << "Touchable not found." << G4endl;
    }

    return;
    
  } else if (command == fpCommandDraw) {

    G4PhysicalVolumeModel::TouchableProperties properties =
    G4TouchableUtils::FindTouchableProperties(fCurrentTouchableProperties.fTouchablePath);
    if (properties.fpTouchablePV) {
      // To handle paramaterisations we have to set the copy number
      properties.fpTouchablePV->SetCopyNo(properties.fCopyNo);
      G4PhysicalVolumeModel* pvModel = new G4PhysicalVolumeModel
      (properties.fpTouchablePV,
       G4PhysicalVolumeModel::UNLIMITED,
       properties.fTouchableGlobalTransform,
       nullptr, // Modelling parameters (not used)
       true, // use full extent (prevents calculating own extent, which crashes)
       properties.fTouchableBaseFullPVPath);

      UImanager->ApplyCommand("/vis/scene/create");
      currentScene = fpVisManager->GetCurrentScene();  // New current scene
      G4bool successful = currentScene->AddRunDurationModel(pvModel,warn);
      UImanager->ApplyCommand("/vis/sceneHandler/attach");

      if (successful) {
        if (fpCommandDraw->GetNewBoolValue(newValue)) {
          const auto& extent = pvModel->GetExtent();
          const G4double halfX = (extent.GetXmax()-extent.GetXmin())/2.;
          const G4double halfY = (extent.GetYmax()-extent.GetYmin())/2.;
          const G4double halfZ = (extent.GetZmax()-extent.GetZmin())/2.;
          G4Box extentBox("extent",halfX,halfY,halfZ);
          G4VisAttributes extentVA;
          extentVA.SetForceWireframe();
          fpVisManager->Draw(extentBox,extentVA,G4Translate3D(extent.GetExtentCentre()));
        }
        if (verbosity >= G4VisManager::confirmations) {
          G4cout << "\"" << properties.fpTouchablePV->GetName()
          << "\", copy no. " << properties.fCopyNo << " drawn";
          if (fpCommandDraw->GetNewBoolValue(newValue)) {
            G4cout << " with extent box";
          }
          G4cout << '.' << G4endl;
        }
      } else {
        G4VisCommandsSceneAddUnsuccessful(verbosity);
      }
    } else {
      G4warn << "Touchable not found." << G4endl;
    }
    return;

  } else if (command == fpCommandDump) {

    G4PhysicalVolumeModel::TouchableProperties properties =
    G4TouchableUtils::FindTouchableProperties(fCurrentTouchableProperties.fTouchablePath);
    if (properties.fpTouchablePV) {
      // To handle paramaterisations we have to set the copy number
      properties.fpTouchablePV->SetCopyNo(properties.fCopyNo);
      G4PhysicalVolumeModel tempPVModel
      (properties.fpTouchablePV,
       G4PhysicalVolumeModel::UNLIMITED,
       properties.fTouchableGlobalTransform,
       nullptr, // Modelling parameters (not used)
       true, // use full extent (prevents calculating own extent, which crashes)
       properties.fTouchableBaseFullPVPath);
      const std::map<G4String,G4AttDef>* attDefs = tempPVModel.GetAttDefs();
      std::vector<G4AttValue>* attValues = tempPVModel.CreateCurrentAttValues();
      G4cout << G4AttCheck(attValues,attDefs);
      delete attValues;
      const auto lv = properties.fpTouchablePV->GetLogicalVolume();
      const auto polyhedron = lv->GetSolid()->GetPolyhedron();
      polyhedron->SetVisAttributes(lv->GetVisAttributes());
      G4cout << "\nLocal polyhedron coordinates:\n" << *polyhedron;
      const G4Transform3D& transform = tempPVModel.GetCurrentTransform();
      polyhedron->Transform(transform);
      G4cout << "\nGlobal polyhedron coordinates:\n" << *polyhedron;
    } else {
      G4warn << "Touchable not found." << G4endl;
    }
    return;

  } else if (command == fpCommandExtentForField) {

    G4PhysicalVolumeModel::TouchableProperties properties =
    G4TouchableUtils::FindTouchableProperties(fCurrentTouchableProperties.fTouchablePath);
    if (properties.fpTouchablePV) {
      G4VisExtent extent
      = properties.fpTouchablePV->GetLogicalVolume()->GetSolid()->GetExtent();
      extent.Transform(properties.fTouchableGlobalTransform);
      fCurrentExtentForField = extent;
      fCurrrentPVFindingsForField.clear();
      if (verbosity >= G4VisManager::confirmations) {
        G4cout << "Extent for field set to " << extent
        << "\nVolume for field has been cleared."
        << G4endl;
      }
      if (fpCommandExtentForField->GetNewBoolValue(newValue)) {
        DrawExtent(extent);
      }
    } else {
      G4warn << "Touchable not found." << G4endl;
    }
    return;

  } else if (command == fpCommandFindPath) {

    G4String pvName;
    G4int copyNo;
    std::istringstream iss(newValue);
    iss >> pvName >> copyNo;
    std::vector<G4PhysicalVolumesSearchScene::Findings> findingsVector;
    std::vector<G4VPhysicalVolume*>::iterator iterWorld =
    transportationManager->GetWorldsIterator();
    for (size_t i = 0; i < nWorlds; ++i, ++iterWorld) {
      G4PhysicalVolumeModel searchModel (*iterWorld);  // Unlimited depth.
      G4ModelingParameters mp;  // Default - no culling.
      searchModel.SetModelingParameters (&mp);
      // Find all instances at any position in the tree
      G4PhysicalVolumesSearchScene searchScene (&searchModel, pvName, copyNo);
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
      G4warn << pvName;
      if (copyNo >= 0) G4warn << ':' << copyNo;
      G4warn << " not found" << G4endl;
    }

  } else if (command == fpCommandLocalAxes) {

    const auto& transform = fCurrentTouchableProperties.fTouchableGlobalTransform;
    const auto& extent = fCurrentTouchableProperties.fpTouchablePV->GetLogicalVolume()->GetSolid()->GetExtent();
    const G4double lengthMax = extent.GetExtentRadius()/2.;
    const G4double intLog10LengthMax = std::floor(std::log10(lengthMax));
    G4double length = std::pow(10,intLog10LengthMax);
    if (5.*length < lengthMax) length *= 5.;
    else if (2.*length < lengthMax) length *= 2.;
    G4AxesModel axesModel(0.,0.,0.,length,transform);
    axesModel.SetGlobalTag("LocalAxesModel");
    axesModel.DescribeYourselfTo(*fpVisManager->GetCurrentSceneHandler());

  } else if (command == fpCommandShowExtent) {

    G4PhysicalVolumeModel::TouchableProperties properties =
    G4TouchableUtils::FindTouchableProperties(fCurrentTouchableProperties.fTouchablePath);
    if (properties.fpTouchablePV) {
      G4VisExtent extent
      = properties.fpTouchablePV->GetLogicalVolume()->GetSolid()->GetExtent();
      extent.Transform(properties.fTouchableGlobalTransform);
      G4cout << extent << G4endl;
      if (fpCommandShowExtent->GetNewBoolValue(newValue)) DrawExtent(extent);
    } else {
      G4warn << "Touchable not found." << G4endl;
    }
    return;

  } else if (command == fpCommandVolumeForField) {

    G4PhysicalVolumeModel::TouchableProperties properties =
    G4TouchableUtils::FindTouchableProperties(fCurrentTouchableProperties.fTouchablePath);
    if (properties.fpTouchablePV) {
      G4VisExtent extent
      = properties.fpTouchablePV->GetLogicalVolume()->GetSolid()->GetExtent();
      extent.Transform(properties.fTouchableGlobalTransform);
      fCurrentExtentForField = extent;
      fCurrrentPVFindingsForField.clear();
      fCurrrentPVFindingsForField.push_back
      (G4PhysicalVolumesSearchScene::Findings(properties));
      if (verbosity >= G4VisManager::confirmations) {
        G4cout
        << "Volume for field set to " << properties.fpTouchablePV->GetName()
        << ':' << properties.fCopyNo
        << " at " << properties.fTouchableBaseFullPVPath
        << G4endl;
      }
      if (fpCommandVolumeForField->GetNewBoolValue(newValue)) {
        DrawExtent(extent);
      }
    } else {
      G4warn << "Touchable not found." << G4endl;
    }
    return;

  } else {

    if (verbosity >= G4VisManager::errors) {
      G4warn <<
      "ERROR: G4VisCommandsTouchable::SetNewValue: unrecognised command."
      << G4endl;
    }
    return;
  }
}

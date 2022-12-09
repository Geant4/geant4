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

// Base class for visualization commands - John Allison  9th August 1998
// It is really a messenger - we have one command per messenger.

#include "G4VVisCommand.hh"

#include "G4UIcommand.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include <sstream>
#include <cctype>

#include "G4PhysicalVolumeModel.hh"
#include "G4LogicalVolume.hh"

#define G4warn G4cout

G4int           G4VVisCommand::fCurrentArrow3DLineSegmentsPerCircle = 6;
G4Colour        G4VVisCommand::fCurrentColour = G4Colour::White();
G4Colour        G4VVisCommand::fCurrentTextColour = G4Colour::Blue();
G4Text::Layout  G4VVisCommand::fCurrentTextLayout = G4Text::left;
G4double        G4VVisCommand::fCurrentTextSize = 12.;  // pixels
G4double        G4VVisCommand::fCurrentLineWidth = 1.;  // pixels
// Not yet used: G4VisAttributes::LineStyle G4VVisCommand::fCurrentLineStyle = G4VisAttributes::unbroken;
// Not yet used: G4VMarker::FillStyle       G4VVisCommand::fCurrentFillStyle = G4VMarker::filled;
// Not yet used: G4VMarker::SizeType        G4VVisCommand::fCurrentSizeType = G4VMarker::screen;
G4PhysicalVolumeModel::TouchableProperties          G4VVisCommand::fCurrentTouchableProperties;
G4VisExtent                                         G4VVisCommand::fCurrentExtentForField;
std::vector<G4PhysicalVolumesSearchScene::Findings> G4VVisCommand::fCurrrentPVFindingsForField;
G4bool G4VVisCommand::fThereWasAViewer = false;
G4ViewParameters G4VVisCommand::fExistingVP;

G4VVisCommand::G4VVisCommand () {}

G4VVisCommand::~G4VVisCommand () {}

G4VisManager* G4VVisCommand::fpVisManager = nullptr;

G4VisManager* G4VVisCommand::GetVisManager ()
{
  return fpVisManager;
}

void G4VVisCommand::SetVisManager (G4VisManager* pVisManager)
{
  fpVisManager = pVisManager;
}

const G4Colour& G4VVisCommand::GetCurrentTextColour()
{
  return fCurrentTextColour;
}

G4String G4VVisCommand::ConvertToString
(G4double x, G4double y, const char * unitName)
{
  G4double uv = G4UIcommand::ValueOf(unitName);
  
  std::ostringstream oss;
  oss << x/uv << " " << y/uv << " " << unitName;
  return oss.str();
}

G4bool G4VVisCommand::ConvertToDoublePair(const G4String& paramString,
					G4double& xval,
					G4double& yval)
{
  G4double x, y;
  G4String unit;
  
  std::istringstream is(paramString);
  is >> x >> y >> unit;

  if (G4UnitDefinition::IsUnitDefined(unit)) {
    xval = x*G4UIcommand::ValueOf(unit);
    yval = y*G4UIcommand::ValueOf(unit);
  } else {
    G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
    if (verbosity >= G4VisManager::errors) {
      G4warn << "ERROR: Unrecognised unit" << G4endl;
    }
    return false;
  }

  return true;
}

const G4String& G4VVisCommand::ConvertToColourGuidance()
{
  static G4String guidance
  ("Accepts (a) RGB triplet. e.g., \".3 .4 .5\", or"
   "\n (b) string such as \"white\", \"black\", \"grey\", \"red\"...or"
   "\n (c) an additional number for opacity, e.g., \".3 .4 .5 .6\""
   "\n     or \"grey ! ! .6\" (note \"!\"'s for unused parameters).");
  return guidance;
}

void G4VVisCommand::ConvertToColour
(G4Colour& colour,
 const G4String& redOrString, G4double green, G4double blue, G4double opacity)
{
  // Note: colour is supplied by the caller and some or all of its components
  // may act as default.
  //
  // Note: redOrString is either a number or string.  If a string it must be
  // one of the recognised colours.
  //
  // Thus the arguments can be, for example:
  // (colour,"red",...,...,0.5): will give the colour red with opacity 0.5 (the
  // third and fourth arguments are ignored), or
  // (1.,0.,0.,0.5): this also will be red with opacity 0.5.

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  const std::size_t iPos0 = 0;
  if (std::isalpha(redOrString[iPos0])) {

    // redOrString is probably alphabetic characters defining the colour
    if (!G4Colour::GetColour(redOrString, colour)) {
      // Not a recognised string
      if (verbosity >= G4VisManager::warnings) {
        G4warn << "WARNING: Colour \"" << redOrString
        << "\" not found.  Defaulting to " << colour
        << G4endl;
      }
      return;
    } else {
      // It was a recognised string.  Now add opacity.
      colour.SetAlpha(opacity);
      return;
    }

  } else {

    // redOrString is probably numeric defining the red component
    std::istringstream iss(redOrString);
    G4double red;
    iss >> red;
    if (iss.fail()) {
      if (verbosity >= G4VisManager::warnings) {
        G4warn << "WARNING: String \"" << redOrString
        << "\" cannot be parsed.  Defaulting to " << colour
        << G4endl;
      }
      return;
    } else {
      colour = G4Colour(red,green,blue,opacity);
      return;
    }
    
  }
}

G4bool G4VVisCommand::ProvideValueOfUnit
(const G4String& where,
 const G4String& unit,
 const G4String& category,
 G4double& value)
{
  // Return false if there's a problem

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4bool success = true;
  if (!G4UnitDefinition::IsUnitDefined(unit)) {
    if (verbosity >= G4VisManager::warnings) {
      G4warn << where
      << "\n  Unit \"" << unit << "\" not defined"
      << G4endl;
    }
    success = false;
  } else if (G4UnitDefinition::GetCategory(unit) != category) {
    if (verbosity >= G4VisManager::warnings) {
      G4warn << where
      << "\n  Unit \"" << unit << "\" not a unit of " << category;
      if (category == "Volumic Mass") G4warn << " (density)";
      G4warn << G4endl;
    }
    success = false;
  } else {
    value = G4UnitDefinition::GetValueOf(unit);
  }
  return success;
}

void G4VVisCommand::CheckSceneAndNotifyHandlers(G4Scene* pScene)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  if (!pScene) {
    if (verbosity >= G4VisManager::warnings) {
      G4warn << "WARNING: Scene pointer is null."
      << G4endl;
    }
    return;
  }

  G4VSceneHandler* pSceneHandler = fpVisManager -> GetCurrentSceneHandler();
  if (!pSceneHandler) {
    if (verbosity >= G4VisManager::warnings) {
      G4warn << "WARNING: Scene handler not found." << G4endl;
    }
    return;
  }

  // Scene has changed.  If it is the scene of the currrent scene handler
  // refresh viewers of all scene handlers using this scene. If not, it may be
  // a scene that the user is building up before attaching to a scene handler,
  // so do nothing.
  if (pScene == pSceneHandler->GetScene()) {
    G4UImanager::GetUIpointer () -> ApplyCommand ("/vis/scene/notifyHandlers");
  }

}

G4bool G4VVisCommand::CheckView ()
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();

  if (!viewer) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<
  "ERROR: No current viewer - \"/vis/viewer/list\" to see possibilities."
       << G4endl;
    }
    return false;
  }
  
  return true;
}

void G4VVisCommand::G4VisCommandsSceneAddUnsuccessful
(G4VisManager::Verbosity verbosity) {
  // Some frequently used error printing...
  if (verbosity >= G4VisManager::warnings) {
    G4warn <<
    "WARNING: For some reason, possibly mentioned above, it has not been"
    "\n  possible to add to the scene."
    << G4endl;
  }
}

void G4VVisCommand::SetViewParameters
(G4VViewer* viewer, const G4ViewParameters& viewParams) {
  viewer->SetViewParameters(viewParams);
  RefreshIfRequired(viewer);
}

void G4VVisCommand::RefreshIfRequired(G4VViewer* viewer) {
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4VSceneHandler* sceneHandler = viewer->GetSceneHandler();
  const G4ViewParameters& viewParams = viewer->GetViewParameters();
  if (sceneHandler && sceneHandler->GetScene()) {
    if (viewParams.IsAutoRefresh()) {
      G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/refresh");
    }
    else {
      if (verbosity >= G4VisManager::warnings) {
        G4warn << "Issue /vis/viewer/refresh or flush to see effect."
        << G4endl;
      }
    }
  }
}

void G4VVisCommand::InterpolateViews
(G4VViewer* currentViewer,
 std::vector<G4ViewParameters> viewVector,
 const G4int nInterpolationPoints,
 const G4int waitTimePerPointmilliseconds,
 const G4String exportString)
{
  const G4int safety = (G4int)viewVector.size()*nInterpolationPoints;
  G4int safetyCount = 0;
  do {
    G4ViewParameters* vp =
    G4ViewParameters::CatmullRomCubicSplineInterpolation(viewVector,nInterpolationPoints);
    if (!vp) break;  // Finished.
    currentViewer->SetViewParameters(*vp);
    currentViewer->RefreshView();
    if (exportString == "export" &&
        G4StrUtil::contains(currentViewer->GetName(), "OpenGL")) {
      G4UImanager::GetUIpointer()->ApplyCommand("/vis/ogl/export");
    }
    currentViewer->ShowView();
    if (waitTimePerPointmilliseconds > 0)
      std::this_thread::sleep_for(std::chrono::milliseconds(waitTimePerPointmilliseconds));
  } while (safetyCount++ < safety);  // Loop checking, 16.02.2016, J.Allison
}

void G4VVisCommand::InterpolateToNewView
(G4VViewer* currentViewer,
 const G4ViewParameters& oldVP,
 const G4ViewParameters& newVP,
 const G4int nInterpolationPoints,
 const G4int waitTimePerPointmilliseconds,
 const G4String exportString)
{
  std::vector<G4ViewParameters> viewVector;
  viewVector.push_back(oldVP);
  viewVector.push_back(oldVP);
  viewVector.push_back(newVP);
  viewVector.push_back(newVP);

  InterpolateViews
  (currentViewer,
   viewVector,
   nInterpolationPoints,
   waitTimePerPointmilliseconds,
   exportString);
}

void G4VVisCommand::Twinkle
// Twinkles the touchables in paths
// /vis/viewer/centreOn to see its effect
 (G4VViewer* currentViewer,
  const G4ViewParameters& baseVP,
  const std::vector<std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>>& paths)
{
  // Copy view parameters to temporary variables ready for adding VisAttributes Modifiers (VAMs)
  auto loVP = baseVP;  // For black and solid VAMs
  auto hiVP = baseVP;  // For white and solid VAMs

  // Modify them with vis attribute modifiers (VAMs)
  for (const auto& path: paths) {
    const auto& touchable = path.back().GetPhysicalVolume();
    auto loVisAtts
    = *(currentViewer->GetApplicableVisAttributes
        (touchable->GetLogicalVolume()->GetVisAttributes()));
    auto hiVisAtts = loVisAtts;
    loVisAtts.SetColour(G4Colour::Black());
    loVisAtts.SetForceSolid();
    hiVisAtts.SetColour(G4Colour::White());
    hiVisAtts.SetForceSolid();
    auto pvNameCopyNoPath
    = G4PhysicalVolumeModel::GetPVNameCopyNoPath(path);

    auto loVAMColour = G4ModelingParameters::VisAttributesModifier
    (loVisAtts, G4ModelingParameters::VASColour, pvNameCopyNoPath);
    loVP.AddVisAttributesModifier(loVAMColour);
    auto loVAMStyle = G4ModelingParameters::VisAttributesModifier
    (loVisAtts, G4ModelingParameters::VASForceSolid, pvNameCopyNoPath);
    loVP.AddVisAttributesModifier(loVAMStyle);

    auto hiVAMColour = G4ModelingParameters::VisAttributesModifier
    (hiVisAtts, G4ModelingParameters::VASColour, pvNameCopyNoPath);
    hiVP.AddVisAttributesModifier(hiVAMColour);
    auto hiVAMStyle = G4ModelingParameters::VisAttributesModifier
    (hiVisAtts, G4ModelingParameters::VASForceSolid, pvNameCopyNoPath);
    hiVP.AddVisAttributesModifier(hiVAMStyle);
  }

  // Twinkle
  std::vector<G4ViewParameters> viewVector;
  viewVector.push_back(loVP);
  viewVector.push_back(hiVP);
  viewVector.push_back(loVP);
  viewVector.push_back(hiVP);
  viewVector.push_back(loVP);
  viewVector.push_back(hiVP);
  viewVector.push_back(loVP);
  viewVector.push_back(hiVP);
  viewVector.push_back(loVP);
  viewVector.push_back(hiVP);
  // Just 5 interpolation points for a reasonable twinkle rate
  InterpolateViews(currentViewer,viewVector,5);
}

void G4VVisCommand::CopyGuidanceFrom
(const G4UIcommand* fromCmd, G4UIcommand* toCmd, G4int startLine)
{
  if (fromCmd && toCmd) {
    const G4int nGuideEntries = (G4int)fromCmd->GetGuidanceEntries();
    for (G4int i = startLine; i < nGuideEntries; ++i) {
      const G4String& guidance = fromCmd->GetGuidanceLine(i);
      toCmd->SetGuidance(guidance);
    }
  }
}

void G4VVisCommand::CopyParametersFrom
(const G4UIcommand* fromCmd, G4UIcommand* toCmd)
{
  if (fromCmd && toCmd) {
    const G4int nParEntries = (G4int)fromCmd->GetParameterEntries();
    for (G4int i = 0; i < nParEntries; ++i) {
      G4UIparameter* parameter = new G4UIparameter(*(fromCmd->GetParameter(i)));
      toCmd->SetParameter(parameter);
    }
  }
}

void G4VVisCommand::DrawExtent(const G4VisExtent& extent)
{
  if (fpVisManager) {
    const G4double halfX = (extent.GetXmax() - extent.GetXmin()) / 2.;
    const G4double halfY = (extent.GetYmax() - extent.GetYmin()) / 2.;
    const G4double halfZ = (extent.GetZmax() - extent.GetZmin()) / 2.;
    if (halfX > 0. && halfY > 0. && halfZ > 0.) {
      const G4Box box("vis_extent",halfX,halfY,halfZ);
      const G4VisAttributes visAtts(G4Color::Red());
      const G4Point3D& centre = extent.GetExtentCenter();
      fpVisManager->Draw(box,visAtts,G4Translate3D(centre));
    }
  }
}

void G4VVisCommand::CopyCameraParameters
(G4ViewParameters& target, const G4ViewParameters& from)
{
  // Copy view parameters pertaining only to camera
  target.SetViewpointDirection  (from.GetViewpointDirection());
  target.SetLightpointDirection (from.GetLightpointDirection());
  target.SetLightsMoveWithCamera(from.GetLightsMoveWithCamera());
  target.SetUpVector            (from.GetUpVector());
  target.SetFieldHalfAngle      (from.GetFieldHalfAngle());
  target.SetZoomFactor          (from.GetZoomFactor());
  target.SetScaleFactor         (from.GetScaleFactor());
  target.SetCurrentTargetPoint  (from.GetCurrentTargetPoint());
  target.SetDolly               (from.GetDolly());
}

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
// $Id: G4VVisCommand.cc 104163 2017-05-15 06:52:42Z gcosmo $

// Base class for visualization commands - John Allison  9th August 1998
// It is really a messenger - we have one command per messenger.

#include "G4VVisCommand.hh"

#include "G4UIcommand.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include <sstream>
#include <cctype>

G4int G4VVisCommand::fErrorCode = 0;

G4int           G4VVisCommand::fCurrentArrow3DLineSegmentsPerCircle = 6;
G4Colour        G4VVisCommand::fCurrentColour = G4Colour::White();
G4Colour        G4VVisCommand::fCurrentTextColour = G4Colour::Blue();
G4Text::Layout  G4VVisCommand::fCurrentTextLayout = G4Text::left;
G4double        G4VVisCommand::fCurrentTextSize = 12.;  // pixels
G4double        G4VVisCommand::fCurrentLineWidth = 1.;  // pixels
// Not yet used: G4VisAttributes::LineStyle G4VVisCommand::fCurrentLineStyle = G4VisAttributes::unbroken;
// Not yet used: G4VMarker::FillStyle       G4VVisCommand::fCurrentFillStyle = G4VMarker::filled;
// Not yet used: G4VMarker::SizeType        G4VVisCommand::fCurrentSizeType = G4VMarker::screen;
G4ModelingParameters::PVNameCopyNoPath G4VVisCommand::fCurrentTouchablePath;

G4VVisCommand::G4VVisCommand () {}

G4VVisCommand::~G4VVisCommand () {}

G4VisManager* G4VVisCommand::fpVisManager = 0;

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
      G4cout << "ERROR: Unrecognised unit" << G4endl;
    }
    return false;
  }

  return true;
}

const G4String& G4VVisCommand::ConvertToColourGuidance()
{
  static G4String guidance
  ("Accepts (a) RGB triplet. e.g., \".3 .4 .5\", or"
   "\n(b) string such as \"white\", \"black\", \"grey\", \"red\"..."
   "\n(c) an additional number for opacity, e.g., \".3 .4 .5 .6\""
   "\n    or \"grey ! ! .6\" (note \"!\"'s for unused green and blue parameters),"
   "\n    e.g. \"! ! ! 0.\" for a transparent colour.");
  return guidance;
}

void G4VVisCommand::ConvertToColour
(G4Colour& colour,
 const G4String& redOrString, G4double green, G4double blue, G4double opacity)
{
  // Note: colour is supplied by the caller and becomes the default if the
  // remaining parameters cannot be parsed.

  // Note: redOrString is either a number or string.  If a string it must be
  // one of the recognised colours.

  // Thus the arguments can be, for example:
  // (colour,"red",...,...,0.5): will give the colour red with opacity 0.5 (the
  // third and fourth arguments are ignored), or
  // (1.,0.,0.,0.5): this also will be red with opacity 0.5.

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  const size_t iPos0 = 0;
  if (std::isalpha(redOrString[iPos0])) {

    // redOrString is probably alphabetic characters defining the colour
    if (!G4Colour::GetColour(redOrString, colour)) {
      // Not a recognised string
      if (verbosity >= G4VisManager::warnings) {
        G4cout << "WARNING: Colour \"" << redOrString
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
        G4cout << "WARNING: String \"" << redOrString
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

void G4VVisCommand::UpdateVisManagerScene
(const G4String& sceneName) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  const G4SceneList& sceneList = fpVisManager -> GetSceneList ();

  G4int iScene, nScenes = sceneList.size ();
  for (iScene = 0; iScene < nScenes; iScene++) {
    if (sceneList [iScene] -> GetName () == sceneName) break;
  }

  G4Scene* pScene = 0;  // Zero unless scene has been found...
  if (iScene < nScenes) {
    pScene = sceneList [iScene];
  }

  if (!pScene) {
    if (verbosity >= G4VisManager::warnings) {
      G4cout << "WARNING: Scene \"" << sceneName << "\" not found."
	     << G4endl;
    }
    return;
  }

  fpVisManager -> SetCurrentScene (pScene);

  // Scene has changed.  Refresh viewers of all sceneHandlers using
  // this scene...
  G4VViewer* pViewer = fpVisManager -> GetCurrentViewer();
  G4VSceneHandler* sceneHandler = fpVisManager -> GetCurrentSceneHandler();
  if (sceneHandler && sceneHandler -> GetScene ()) {
    if (pViewer) {
      G4UImanager::GetUIpointer () ->
	ApplyCommand ("/vis/scene/notifyHandlers");
    }
  }
}

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

#ifndef G4VVISCOMMAND_HH
#define G4VVISCOMMAND_HH

#include "G4VisManager.hh"
#include "G4UImessenger.hh"
#include "G4ThreeVector.hh"
#include "G4Text.hh"
#include "G4VisAttributes.hh"
#include "G4VMarker.hh"
#include "G4ModelingParameters.hh"
#include "G4PhysicalVolumesSearchScene.hh"
#include <vector>

class G4UIcommand;
class G4UIcmdWithAString;

class G4VVisCommand: public G4UImessenger
{
public:
  
  // Uses compiler defaults for copy constructor and assignment.
  G4VVisCommand ();
  virtual ~G4VVisCommand ();

  static G4VisManager* GetVisManager ();

  static void SetVisManager (G4VisManager* pVisManager);

  static const G4Colour& GetCurrentTextColour();

protected:

  // Utility functions

  void SetViewParameters(G4VViewer* viewer, const G4ViewParameters& viewParams);

  void RefreshIfRequired(G4VViewer* viewer);

  void InterpolateViews
  (G4VViewer* currentViewer,
   std::vector<G4ViewParameters> viewVector,
   const G4int nInterpolationPoints = 50,
   const G4int waitTimePerPointmilliseconds = 20,
   const G4String exportString = "");

  void InterpolateToNewView
  (G4VViewer* currentViewer,
   const G4ViewParameters& oldVP,
   const G4ViewParameters& newVP,
   const G4int nInterpolationPoints = 50,
   const G4int waitTimePerPointmilliseconds = 20,
   const G4String exportString = "");

  void Twinkle
  // Twinkles the touchables in paths
  // /vis/viewer/centreOn to see its effect
  (G4VViewer* currentViewer,
   const G4ViewParameters& baseVP,
   const std::vector<std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>>& paths);

  // Conversion routines augmenting those in G4UIcommand.

  static G4String ConvertToString(G4double x, G4double y,
				  const char * unitName);

  static G4bool ConvertToDoublePair(const G4String& paramString,
                                    G4double& xval,
                                    G4double& yval);
  // Return false if problem parsing paramString.

  const G4String& ConvertToColourGuidance();
  void ConvertToColour
  (G4Colour& colour,
   const G4String& redOrString,
   G4double green,
   G4double blue,
   G4double opacity);
  // Note: colour is supplied by the caller and becomes the default if the
  // remaining parameters cannot be parsed.
  // Note: redOrString is either a number or string.  If a string it must be
  // one of the recognised colours.
  // Thus the arguments can be, for example:
  // (colour,"red",...,...,0.5): will give the colour red with opacity 0.5 (the
  // third and fourth arguments are ignored), or
  // (1.,0.,0.,0.5): this also will be red with opacity 0.5.

  G4bool ProvideValueOfUnit
  (const G4String& where,
   const G4String& unit,
   const G4String& category,
   G4double& value);
  // Return false if there's a problem

  void CopyCameraParameters
  (G4ViewParameters& target, const G4ViewParameters& from);
  // Copy view parameters pertaining only to camera

  // Other utilities

  void CheckSceneAndNotifyHandlers (G4Scene* = nullptr);

  G4bool CheckView();  // False if not valid

  void G4VisCommandsSceneAddUnsuccessful(G4VisManager::Verbosity verbosity);

  void CopyGuidanceFrom
  (const G4UIcommand* fromCmd, G4UIcommand* toCmd, G4int startLine = 0);

  void CopyParametersFrom
  (const G4UIcommand* fromCmd, G4UIcommand* toCmd);

  void DrawExtent(const G4VisExtent&);

  // Data members

  static G4VisManager* fpVisManager;

  // Current quantities for use in appropriate commands
  static G4int fCurrentArrow3DLineSegmentsPerCircle;
  static G4Colour                   fCurrentColour;
  static G4double                   fCurrentLineWidth;
  //static G4VisAttributes::LineStyle fCurrentLineStyle;  Not yet used.
  //static G4VMarker::FillStyle       fCurrentFillStyle;  Not yet used.
  //static G4VMarker::SizeType        fCurrentSizeType;  Not yet used.
  static G4Colour                   fCurrentTextColour;
  static G4Text::Layout             fCurrentTextLayout;
  static G4double                   fCurrentTextSize;
  static G4PhysicalVolumeModel::TouchableProperties fCurrentTouchableProperties;
  static G4VisExtent                fCurrentExtentForField;
  static std::vector<G4PhysicalVolumesSearchScene::Findings> fCurrrentPVFindingsForField;

  // When we create a new viewer we would like to use the view parameters of
  // the existing viewer if there was one. This has to be checked at the
  // creation of a new viewer and *also* at the creation of a new scene
  // handler.
  static G4bool fThereWasAViewer;  // True if there was a viewer
  static G4ViewParameters fExistingVP;  // Its view parameters
};

#endif

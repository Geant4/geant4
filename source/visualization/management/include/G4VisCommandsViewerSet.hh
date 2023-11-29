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

// /vis/viewer/set commands - John Allison  16th May 2000

#ifndef G4VISCOMMANDSVIEWERSET_HH
#define G4VISCOMMANDSVIEWERSET_HH

#include "G4VisCommandsViewer.hh"

#include <vector>

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWith3VectorAndUnit;

class G4VisCommandsViewerSet: public G4VVisCommand {
public:
  G4VisCommandsViewerSet ();
  virtual ~G4VisCommandsViewerSet ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandsViewerSet (const G4VisCommandsViewerSet&);
  G4VisCommandsViewerSet& operator = (const G4VisCommandsViewerSet&);
  G4UIcmdWithAString*   fpCommandAll;
  G4UIcmdWithABool*     fpCommandAutoRefresh;
  G4UIcmdWithABool*     fpCommandAuxEdge;
  G4UIcommand*          fpCommandBackground;
  G4UIcommand*          fpCommandCulling;
  G4UIcmdWithAString*   fpCommandCutawayMode;
  G4UIcommand*          fpCommandDefaultColour;
  G4UIcommand*          fpCommandDefaultTextColour;
  G4UIcmdWithABool*     fpCommandEdge;
  G4UIcommand*          fpCommandExplodeFactor;
  G4UIcmdWithADouble*   fpCommandGlobalMarkerScale;
  G4UIcmdWithADouble*   fpCommandGlobalLineWidthScale;
  G4UIcmdWithABool*     fpCommandHiddenEdge;
  G4UIcmdWithABool*     fpCommandHiddenMarker;
  G4UIcmdWithAString*   fpCommandLightsMove;
  G4UIcommand*          fpCommandLightsThetaPhi;
  G4UIcommand*          fpCommandLightsVector;
  G4ThreeVector         fLightsVector;
  G4UIcmdWithAnInteger* fpCommandLineSegments;
  G4UIcmdWithoutParameter* fpCommandLineWidth;
  G4UIcmdWithAnInteger* fpCommandNumberOfCloudPoints;
  G4UIcmdWithABool*     fpCommandPicking;
  G4UIcommand*          fpCommandProjection;
  G4UIcmdWithAString*   fpCommandRotationStyle;
  G4UIcommand*          fpCommandSectionPlane;
  G4UIcmdWithABool*     fpCommandSpecialMeshRendering;
  G4UIcmdWithAString*   fpCommandSpecialMeshRenderingOption;
  G4UIcommand*          fpCommandSpecialMeshVolumes;
  G4UIcmdWithAString*   fpCommandStyle;
  G4UIcmdWith3VectorAndUnit* fpCommandTargetPoint;
  G4UIcommand*          fpCommandUpThetaPhi;
  G4UIcommand*          fpCommandUpVector;
  G4ThreeVector         fUpVector;
  G4UIcommand*          fpCommandViewpointThetaPhi;
  G4UIcommand*          fpCommandViewpointVector;
  G4ThreeVector         fViewpointVector;
  G4UIdirectory*        fpTimeWindowDirectory;
  G4UIcommand*          fpCommandTimeWindowDisplayHeadTime;
  G4UIcommand*          fpCommandTimeWindowDisplayLightFront;
  G4UIcommand*          fpCommandTimeWindowEndTime;
  G4UIcmdWithADouble*   fpCommandTimeWindowFadeFactor;
  G4UIcommand*          fpCommandTimeWindowStartTime;
};

#endif

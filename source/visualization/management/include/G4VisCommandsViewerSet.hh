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
// $Id: G4VisCommandsViewerSet.hh,v 1.11 2002-11-27 12:33:29 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// /vis/viewer/set commands - John Allison  16th May 2000

#ifndef G4VISCOMMANDSVIEWERSET_HH
#define G4VISCOMMANDSVIEWERSET_HH

#include "G4VisCommandsViewer.hh"

#include "g4std/vector"

class G4UIcommand;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;

class G4VisCommandsViewerSet: public G4VVisCommandViewer {
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
  G4UIcommand*          fpCommandCulling;
  G4UIcmdWithABool*     fpCommandEdge;
  G4UIcmdWithADouble*   fpCommandGlobalMarkerScale;
  G4UIcmdWithABool*     fpCommandHiddenEdge;
  G4UIcmdWithABool*     fpCommandHiddenMarker;
  G4UIcmdWithAnInteger* fpCommandLineSegments;
  G4UIcmdWithAString*   fpCommandLightsMove;
  G4UIcommand*          fpCommandLightsThetaPhi;
  G4UIcommand*          fpCommandLightsVector;
  G4ThreeVector         fLightsVector;
  G4UIcommand*          fpCommandProjection;
  G4UIcommand*          fpCommandSectionPlane;
  G4UIcmdWithAString*   fpCommandStyle;
  G4UIcommand*          fpCommandUpThetaPhi;
  G4UIcommand*          fpCommandUpVector;
  G4ThreeVector         fUpVector;
  G4UIcommand*          fpCommandViewpointThetaPhi;
  G4UIcommand*          fpCommandViewpointVector;
  G4ThreeVector         fViewpointVector;
};

#endif

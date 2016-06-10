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
// $Id: G4VisCommandsTouchableSet.hh 66373 2012-12-18 09:41:34Z gcosmo $

// /vis/touchable/set commands - John Allison  8th October 2012

#ifndef G4VISCOMMANDSTOUCHABLESET_HH
#define G4VISCOMMANDSTOUCHABLESET_HH

#include "G4VisCommandsViewer.hh"

class G4UIcommand;
class G4UIcmdWithABool;
class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;

class G4VisCommandsTouchableSet: public G4VVisCommandViewer {
public:
  G4VisCommandsTouchableSet ();
  virtual ~G4VisCommandsTouchableSet ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandsTouchableSet (const G4VisCommandsTouchableSet&);
  G4VisCommandsTouchableSet& operator = (const G4VisCommandsTouchableSet&);
  G4UIcommand*          fpCommandSetColour;
  G4UIcmdWithABool*     fpCommandSetDaughtersInvisible;
  G4UIcmdWithABool*     fpCommandSetForceAuxEdgeVisible;
  G4UIcmdWithAnInteger* fpCommandSetLineSegmentsPerCircle;
  G4UIcmdWithABool*     fpCommandSetForceSolid;
  G4UIcmdWithABool*     fpCommandSetForceWireframe;
  G4UIcmdWithAString*   fpCommandSetLineStyle;
  G4UIcmdWithADouble*   fpCommandSetLineWidth;
  G4UIcmdWithABool*     fpCommandSetVisibility;
};

#endif

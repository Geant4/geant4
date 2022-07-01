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

// /vis/touchable/set commands - John Allison  14th May 2014

#ifndef G4VISCOMMANDSTOUCHABLE_HH
#define G4VISCOMMANDSTOUCHABLE_HH

#include "G4VisCommands.hh"

class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;

class G4VisCommandsTouchable: public G4VVisCommand {
public:
  G4VisCommandsTouchable ();
  virtual ~G4VisCommandsTouchable ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandsTouchable (const G4VisCommandsTouchable&);
  G4VisCommandsTouchable& operator = (const G4VisCommandsTouchable&);
  G4UIcmdWithoutParameter* fpCommandCentreOn;
  G4UIcmdWithoutParameter* fpCommandCentreAndZoomInOn;
  G4UIcmdWithABool*        fpCommandDraw;
  G4UIcmdWithoutParameter* fpCommandDump;
  G4UIcmdWithABool*        fpCommandExtentForField;
  G4UIcommand*             fpCommandFindPath;
  G4UIcmdWithoutParameter* fpCommandLocalAxes;
  G4UIcmdWithABool*        fpCommandShowExtent;
  G4UIcmdWithABool*        fpCommandVolumeForField;
};

#endif

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
// $Id: G4VisCommands.hh,v 1.3 2001-08-05 02:29:02 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// /vis/ top level commands - John Allison  5th February 2001

#ifndef G4VISCOMMANDS_HH
#define G4VISCOMMANDS_HH

#include "G4VVisCommand.hh"

class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

class G4VisCommandEnable: public G4VVisCommand {
public:
  G4VisCommandEnable ();
  virtual ~G4VisCommandEnable ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandEnable (const G4VisCommandEnable&);
  G4VisCommandEnable& operator = (const G4VisCommandEnable&);
  G4UIcmdWithABool* fpCommand;
  G4UIcmdWithoutParameter* fpCommand1;
};

class G4VisCommandVerbose: public G4VVisCommand {
public:
  G4VisCommandVerbose ();
  virtual ~G4VisCommandVerbose ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandVerbose (const G4VisCommandVerbose&);
  G4VisCommandVerbose& operator = (const G4VisCommandVerbose&);
  G4UIcmdWithAString* fpCommand;
};

#endif

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

// /vis/ top level commands - John Allison  5th February 2001

#ifndef G4VISCOMMANDS_HH
#define G4VISCOMMANDS_HH

#include "G4VVisCommand.hh"

class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

class G4VisCommandAbortReviewKeptEvents: public G4VVisCommand {
public:
  G4VisCommandAbortReviewKeptEvents ();
  virtual ~G4VisCommandAbortReviewKeptEvents ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandAbortReviewKeptEvents (const G4VisCommandAbortReviewKeptEvents&);
  G4VisCommandAbortReviewKeptEvents& operator = (const G4VisCommandAbortReviewKeptEvents&);
  G4UIcmdWithABool* fpCommand;
};

class G4VisCommandAbortReviewPlots: public G4VVisCommand {
public:
  G4VisCommandAbortReviewPlots ();
  virtual ~G4VisCommandAbortReviewPlots ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandAbortReviewPlots (const G4VisCommandAbortReviewPlots&);
  G4VisCommandAbortReviewPlots& operator = (const G4VisCommandAbortReviewKeptEvents&);
  G4UIcmdWithABool* fpCommand;
};

class G4VisCommandDrawOnlyToBeKeptEvents: public G4VVisCommand {
public:
  G4VisCommandDrawOnlyToBeKeptEvents ();
  virtual ~G4VisCommandDrawOnlyToBeKeptEvents ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandDrawOnlyToBeKeptEvents (const G4VisCommandDrawOnlyToBeKeptEvents&);
  G4VisCommandDrawOnlyToBeKeptEvents& operator = (const G4VisCommandDrawOnlyToBeKeptEvents&);
  G4UIcmdWithABool* fpCommand;
};

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

class G4VisCommandInitialize: public G4VVisCommand {
public:
  G4VisCommandInitialize ();
  virtual ~G4VisCommandInitialize ();
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandInitialize (const G4VisCommandInitialize&);
  G4VisCommandInitialize& operator = (const G4VisCommandInitialize&);
  G4UIcmdWithoutParameter* fpCommand;
};

class G4VisCommandList: public G4VVisCommand {
public:
  G4VisCommandList ();
  virtual ~G4VisCommandList ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandList (const G4VisCommandList&);
  G4VisCommandList& operator = (const G4VisCommandList&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandReviewKeptEvents: public G4VVisCommand {
public:
  G4VisCommandReviewKeptEvents ();
  virtual ~G4VisCommandReviewKeptEvents ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandReviewKeptEvents (const G4VisCommandReviewKeptEvents&);
  G4VisCommandReviewKeptEvents& operator = (const G4VisCommandReviewKeptEvents&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandReviewPlots: public G4VVisCommand {
public:
  G4VisCommandReviewPlots ();
  virtual ~G4VisCommandReviewPlots ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandReviewPlots (const G4VisCommandReviewPlots&);
  G4VisCommandReviewPlots& operator = (const G4VisCommandReviewPlots&);
  G4UIcmdWithoutParameter* fpCommand;
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

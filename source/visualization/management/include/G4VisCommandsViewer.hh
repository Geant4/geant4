// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsViewer.hh,v 1.1 1999-01-07 16:15:24 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// /vis/viewer commands - John Allison  25th October 1998

#ifndef G4VISCOMMANDSVIEWER_HH
#define G4VISCOMMANDSVIEWER_HH

#include "G4VVisCommand.hh"

class G4UIcommand;
class G4UIcmdWithAString;

class G4VVisCommandViewer: public G4VVisCommand {
public:
  // Uses compiler defaults for destructor, copy constructor and assignment.
  G4VVisCommandViewer ();
protected:
  G4String ShortName (const G4String &);
  void UpdateCandidateLists ();
  static G4String fViewerNameList;
  // member so that it has long life - static because shared between objects.
};

class G4VisCommandViewerCreate: public G4VVisCommandViewer {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandViewerCreate ();
  ~G4VisCommandViewerCreate ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4String NextName ();
  G4UIcommand* fpCommand;
  G4int fId;
};

class G4VisCommandViewerList: public G4VVisCommandViewer {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandViewerList ();
  ~G4VisCommandViewerList ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4UIcommand* fpCommand;
};

class G4VisCommandViewerRemove: public G4VVisCommandViewer {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandViewerRemove ();
  ~G4VisCommandViewerRemove ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandViewerSelect: public G4VVisCommandViewer {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandViewerSelect ();
  ~G4VisCommandViewerSelect ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4UIcmdWithAString* fpCommand;
};

#endif

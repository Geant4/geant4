// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsViewer.hh,v 1.10 2001-02-03 18:39:47 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// /vis/viewer commands - John Allison  25th October 1998

#ifndef G4VISCOMMANDSVIEWER_HH
#define G4VISCOMMANDSVIEWER_HH

#include "G4VVisCommand.hh"

class G4UIcommand;
class G4UIcmdWithAString;

class G4VVisCommandViewer: public G4VVisCommand {
public:
  G4VVisCommandViewer ();
  virtual ~G4VVisCommandViewer ();
protected:
  void UpdateCandidateLists ();
private:
  G4VVisCommandViewer (const G4VVisCommandViewer&);
  G4VVisCommandViewer& operator = (const G4VVisCommandViewer&);
};

class G4VisCommandViewerCreate: public G4VVisCommandViewer {
public:
  G4VisCommandViewerCreate ();
  virtual ~G4VisCommandViewerCreate ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerCreate (const G4VisCommandViewerCreate&);
  G4VisCommandViewerCreate& operator = (const G4VisCommandViewerCreate&);
  G4String NextName ();
  G4UIcommand* fpCommand;
  G4int fId;
};

class G4VisCommandViewerList: public G4VVisCommandViewer {
public:
  G4VisCommandViewerList ();
  virtual ~G4VisCommandViewerList ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerList (const G4VisCommandViewerList&);
  G4VisCommandViewerList& operator = (const G4VisCommandViewerList&);
  G4UIcommand* fpCommand;
};

class G4VisCommandViewerRefresh: public G4VVisCommandViewer {
public:
  G4VisCommandViewerRefresh ();
  virtual ~G4VisCommandViewerRefresh ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerRefresh (const G4VisCommandViewerRefresh&);
  G4VisCommandViewerRefresh& operator = (const G4VisCommandViewerRefresh&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandViewerRemove: public G4VVisCommandViewer {
public:
  G4VisCommandViewerRemove ();
  virtual ~G4VisCommandViewerRemove ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerRemove (const G4VisCommandViewerRemove&);
  G4VisCommandViewerRemove& operator = (const G4VisCommandViewerRemove&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandViewerReset: public G4VVisCommandViewer {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandViewerReset ();
  virtual ~G4VisCommandViewerReset ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerReset (const G4VisCommandViewerReset&);
  G4VisCommandViewerReset& operator = (const G4VisCommandViewerReset&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandViewerSelect: public G4VVisCommandViewer {
public:
  G4VisCommandViewerSelect ();
  virtual ~G4VisCommandViewerSelect ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerSelect (const G4VisCommandViewerSelect&);
  G4VisCommandViewerSelect& operator = (const G4VisCommandViewerSelect&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandViewerUpdate: public G4VVisCommandViewer {
public:
  G4VisCommandViewerUpdate ();
  virtual ~G4VisCommandViewerUpdate ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandViewerUpdate (const G4VisCommandViewerUpdate&);
  G4VisCommandViewerUpdate& operator = (const G4VisCommandViewerUpdate&);
  G4UIcmdWithAString* fpCommand;
  G4UIcmdWithAString* fpCommand1;
};

#endif

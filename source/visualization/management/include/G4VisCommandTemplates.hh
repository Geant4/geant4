// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandTemplates.hh,v 1.2 1999-01-09 16:30:53 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Visualization Command Messenger Templates.
// John Allison  7th September 1997

#ifndef G4VISCOMMANDTEMPLATES_HH
#define G4VISCOMMANDTEMPLATES_HH

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIcommand;

/***************************************
// Vis Command Messenger Template - General Purpose.
template <class T>
class G4VisCommandMessenger: public G4UImessenger {
public:
  G4VisCommandMessenger ();
  ~G4VisCommandMessenger ();
  void SetNewValue (G4UIcommand* command, G4String newValues);
private:
  T fCommand;
};
************************************************/

// Command Messenger Template for directory - OLD STYLE!!
template <class tCommandDirectory>
class G4VisCommandDirectoryMessenger: public G4UImessenger {
public:
  G4VisCommandDirectoryMessenger ();
  ~G4VisCommandDirectoryMessenger ();
private:
  tCommandDirectory fCommandDirectory;
  G4UIcommand* fpCommand;
};

// Command Messenger Template for simple commands with no arguments -
// OLD STYLE!!
template <class tSimpleCommand>
class G4VisSimpleCommandMessenger: public G4UImessenger {
public:
  G4VisSimpleCommandMessenger ();
  ~G4VisSimpleCommandMessenger ();
  void SetNewValue (G4UIcommand* command, G4String newValues);
private:
  tSimpleCommand fSimpleCommand;
  G4UIcommand* fpCommand;
};


// Command Messenger Template for Boolean attributes which can
// usefully be implemented as a button on a GUI - OLD STYLE!!
template <class tButtonCommand>
class G4VisButtonCommandMessenger: public G4UImessenger {
public:
  G4VisButtonCommandMessenger ();
  ~G4VisButtonCommandMessenger ();
  void SetNewValue (G4UIcommand* command, G4String newValues);
private:
  tButtonCommand fButtonCommand;
  G4UIcommand* fpCommand;
};

// OLD STYLE!!
G4UIcommand* G4VisButtonCommandMessengerRegister (G4String commandName,
						  G4String guidance,
						  G4UImessenger* messenger);

// OLD STYLE!!
G4int G4VisButtonCommandMessengerInterpret (G4String newValues);

#include "G4VisCommandTemplates.icc"

#endif

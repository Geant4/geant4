// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIXaw.hh,v 1.4 1999-12-15 14:50:46 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4UIXaw_h
#define G4UIXaw_h 

#if defined(G4UI_BUILD_XAW_SESSION) || defined(G4UI_USE_XAW)

#include <X11/Intrinsic.h>

#include "G4VBasicShell.hh"

class G4UIsession;
class G4UImanager;

// Class description :
//
//  G4UIXaw : class to handle an interactive session.
// G4UIXaw is the "Athena" version of G4UIterminal.
// It has been provided at a time where no free Motif
// where available. Now that some exists (lesstif), 
// G4UIXm should be prefered. The look and feel of 
// G4UIXaw is hugly...
//
//  No command completion is available.
//
// Class description - end :

class G4UIXaw : public G4VBasicShell {
public:
  G4UIXaw(int,char**);
  ~G4UIXaw();
  G4UIsession* SessionStart();
  void Prompt(G4String);
  void SessionTerminate();
  void PauseSessionStart(G4String);
  Widget GetDialog();
private:
  void SecondaryLoop(G4String);
  G4bool GetHelpChoice(G4int&);
  void ExitHelp();
private:
  Widget shell,dialog;
  G4bool fHelp;
  G4int fHelpChoice;
  static void Callback(Widget,XtPointer,XtPointer);
};

#endif

#endif


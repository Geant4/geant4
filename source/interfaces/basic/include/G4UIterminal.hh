// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIterminal.hh,v 1.3 1999-11-08 04:11:12 masayasu Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4UIterminal_h
#define G4UIterminal_h 1

#include "G4VBasicShell.hh"
#include "G4UImanager.hh"
#include <fstream.h>

// class description:
//
// This class inherits the class G4UIsession.
// This is the class to use a character-terminal sesion.

class G4UIterminal : public G4VBasicShell
{
      public: // with description
      G4UIterminal();
      ~G4UIterminal();

      G4UIsession * SessionStart();
      // A character-terminal session  "terminalSession" is instantiated.
      // G4cout stream is redirected by default to the constructed instance.
      // Usage:  G4UIsession * terminalSession = new G4UIterminal;
      // "terminalSession" is started.
      // Usage: terminalSession->SessionStart();
      // "terminalSession"  is deleted.
      // Usage: delete terminalSession;
      //
      void PauseSessionStart(G4String);
      G4int ReceiveG4cout(G4String coutString);
      G4int ReceiveG4cerr(G4String cerrString);
      // These methods are implementation of corresponding virtual methods
      // of G4UIsession class.
  private:
      G4UImanager * UI;
      G4String promptCharacter;
      G4bool iExit;
      G4bool iCont;

  private:
      void ExecuteCommand(G4String);
      G4String GetCommand();
      G4bool GetHelpChoice(G4int&);
      void ExitHelp();
};




#endif


// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIterminal.hh,v 1.1 1999-01-07 16:09:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4UIterminal_h
#define G4UIterminal_h 1

#include "G4VBasicShell.hh"
#include "G4UImanager.hh"
#include <fstream.h>

class G4UIterminal : public G4VBasicShell
{
  public:
      G4UIterminal();
      ~G4UIterminal();

      G4UIsession * SessionStart();
      void PauseSessionStart(G4String);
      G4int ReceiveG4cout(G4String coutString);
      G4int ReceiveG4cerr(G4String cerrString);

  private:
      G4UImanager * UI;
      G4String promptCharacter;
      G4bool iExit;
      G4bool iCont;

  private:
      void ExecuteCommand(G4String);
      G4String GetCommand();
      void ChangeDirectoryCommand(G4String);
      void ListDirectory(G4String);
      void TerminalHelp(G4String);
      void ShowCurrent(G4String);
};




#endif


// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIGAG.hh,v 1.3 1999-11-08 04:14:10 masayasu Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4UIGAG.hh
// GAG(Geant4 adaptive GUI) interface class
// 18.Feb.98 Y.Nagamatu and T.Kodama created G4UIGAG from G4UIterminal

// If you want to use GAG (Tcl and JAVA version),
// read a document file README_Momo.html.

#ifndef G4UIGAG_h
#define G4UIGAG_h 1

#include "G4UIsession.hh"
#include "G4UImanager.hh"
#include <fstream.h>
#include <g4std/vector>

enum UImode { terminal_mode , java_mode, tcl_mode };
enum ChangeOfTree  { notChanged=0, added, deleted, addedAndDeleted };

// class description:
//
// This class inherits the class G4UIsession.
// This is the class which is the front-end to GAG (Geant4 Adaptive GUI).
// Its usage is quite similar to G4UIterminal.

class G4UIGAG : public G4UIsession
{
      public: // with description
      G4UIGAG();
      ~G4UIGAG();

      G4UIsession * SessionStart();

      // A GAG session  "gagSession" is instantiated.
      // G4cout stream is redirected by default to the constructed instance.
      // Usage:  G4UIsession * gagSession = new G4UIGAG;
      // "gagSession" is started.
      // Usage: gagSession->SessionStart();
      // "gagSession"  is deleted.
      // Usage: delete gagSession;
      //
      void PauseSessionStart(G4String);
      G4int ReceiveG4cout(G4String);
      G4int ReceiveG4cerr(G4String);         
      // These methods are implementation of the virtual methods of
      // G4UIsession class. 
      //
      void SessionTerminate();
      void Prompt(G4String);
      G4String GetCommand();
      // These methods are not for users.
  private:
      G4String prefix;
      G4UImanager * UI;
      G4String promptCharacter;
      G4bool iExit;
      G4bool iCont;
      UImode uiMode;

  private:
      G4String JVersion;
      G4String TVersion;
      void ExecuteCommand(G4String);
      void ChangeDirectory(G4String);
      void ListDirectory(G4String);
      void TerminalHelp(G4String);
      G4String ModifyPrefix(G4String);
      G4UIcommandTree * FindDirPath(G4String);
      void ShowCurrent(G4String);
      G4String GetFullPath(G4String);

    void SendCommandProperties(G4UIcommandTree *);
    void SendParameterProperties(G4UIcommandTree *);
    void SendAParamProperty(G4UIcommand *);
    void SendATclParamProperty(G4UIcommand *);
    void CodeGenJavaTree(G4UIcommandTree *,int recursiveLevel); 
    void CodeGenJavaParams(G4UIcommandTree *,int recursiveLevel); 
    void CodeGenTclTree(G4UIcommandTree *, int recursiveLevel); 
    void CodeGenTclParams(G4UIcommandTree *, int recursiveLevel); 
    void SendDisableList(G4UIcommandTree *, int recursiveLevel);

    void NotifyStateChange(void);
    void NotifyCommandUpdate(void);
    void NotifyParameterUpdate(G4UIcommand *);

    int CommandUpdated(void);
    void UpdateState(void);
    void UpdateParamVal(void);  // if param is updated,
                               // call NotifyParameterUpdate()

    // --- the following are used by Notify*Update() and *Updated()
    void GetNewTreeStructure( G4UIcommandTree*,int recursiveLevel);
    void GetNewTreeValues( G4UIcommandTree*,int recursiveLevel);
    
    G4std::vector<G4String> previousTreeCommands;
    G4std::vector<G4String> newTreeCommands;
    G4std::vector<G4String> previousTreeParams;
    G4std::vector<G4String> newTreeParams;
    G4std::vector<G4UIcommand*> previousTreePCP;
    G4std::vector<G4UIcommand*> newTreePCP;

};

#endif

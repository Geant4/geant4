// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIGAG.hh,v 1.1 1999-01-07 16:09:30 gunter Exp $
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

enum UImode { terminal_mode , java_mode, tcl_mode };
enum ChangeOfTree  { notChanged=0, added, deleted, addedAndDeleted };

class G4UIGAG : public G4UIsession
{
  public:
      G4UIGAG();
      ~G4UIGAG();

      G4UIsession * SessionStart();
      void PauseSessionStart(G4String);
      void Prompt(G4String);
      G4String GetCommand();
      G4int ReceiveG4cout(G4String);
      G4int ReceiveG4cerr(G4String);         
      void SessionTerminate();

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
    
    RWTValOrderedVector<G4String> previousTreeCommands;
    RWTValOrderedVector<G4String> newTreeCommands;
    RWTValOrderedVector<G4String> previousTreeParams;
    RWTValOrderedVector<G4String> newTreeParams;
    RWTValOrderedVector<G4UIcommand*> previousTreePCP;
    RWTValOrderedVector<G4UIcommand*> newTreePCP;

};

#endif

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
// $Id: G4UIGAG.hh 66892 2013-01-17 10:57:59Z gunter $
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
#include <fstream>
#include <vector>

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
      void PauseSessionStart(const G4String&);
      G4int ReceiveG4cout(const G4String&);
      G4int ReceiveG4cerr(const G4String&);         
      // These methods are implementation of the virtual methods of
      // G4UIsession class. 
      //
      void SessionTerminate();
      void Prompt(const G4String&);
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
      void ExecuteCommand(const G4String&);
      void ChangeDirectory(const G4String&);
      void ListDirectory(const G4String&);
      void TerminalHelp(const G4String&);
     G4String ModifyPrefix(G4String);
      G4UIcommandTree* FindDirPath(const G4String&);
      void ShowCurrent(const G4String&);
      G4String GetFullPath(const G4String&);

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
    
    std::vector<G4String> previousTreeCommands;
    std::vector<G4String> newTreeCommands;
    std::vector<G4String> previousTreeParams;
    std::vector<G4String> newTreeParams;
    std::vector<G4UIcommand*> previousTreePCP;
    std::vector<G4UIcommand*> newTreePCP;

};

#endif

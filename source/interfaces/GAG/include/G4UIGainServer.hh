//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4UIGainServer.hh,v 1.2 2002-06-12 15:27:23 gcosmo Exp $
// $Name: not supported by cvs2svn $
//


#ifndef G4UIGainServer_h
#define G4UIGainServer_h 1

#include <stdio.h>
#include <sys/types.h>
#include <sys/socket.h>
//#include <sys/un.h>
#include <netinet/in.h>
#include <unistd.h>

#include "g4std/iostream"
#include "g4std/fstream"
#include "globals.hh"
#include "G4UImanager.hh"
#include "G4VBasicShell.hh"

#define SOCK_NAME "/tmp/socket"
#define DEFAULT_PORT 40000; 

enum UImode { terminal_mode , java_mode, tcl_mode };
enum ChangeOfTree  { notChanged=0, added, deleted, addedAndDeleted };

// class description:
//
// This class is similar to G4UIGAG.
// While G4UIGAG provides direct connection to GAG via pipe,
// G4UIGainServer provides socket connection to remote Gain client.
// Gain = Geant4 adaptive interface for network
// Its usage is same as G4UIGAG, except that Geant4 application waits
// a connection from a Gain client.

class G4UIGainServer : public  G4VBasicShell {
private:
  G4String prefix;
  G4String promptCharacter;
  G4UImanager* UI;
  UImode uiMode;
  G4String JVersion;
  G4String TVersion;

  // shell
  //  G4VUIshell* shell;

  // program states
  G4bool iExit;
  G4bool iCont;

  // need for socket
  int socketD[3];
  int port;
  struct sockaddr_in saddr;
  struct sockaddr_in caddr;
  int len;
  int ret;
  char buf[1024];

  // --- the following are used by Notify*Update() and *Updated()
  void GetNewTreeStructure( G4UIcommandTree*,int recursiveLevel);
  void GetNewTreeValues( G4UIcommandTree*,int recursiveLevel);
  
  G4std::vector<G4String> previousTreeCommands;
  G4std::vector<G4String> newTreeCommands;
  G4std::vector<G4String> previousTreeParams;
  G4std::vector<G4String> newTreeParams;
  G4std::vector<G4UIcommand*> previousTreePCP;
  G4std::vector<G4UIcommand*> newTreePCP;

public:
  // These methods are implementation of corresponding virtual methods
  // of G4UIsession class.

  // A GainServer session  "gainSession" is instantiated.
  // G4cout stream is redirected by default to the constructed instance.
  // Usage:  G4UIsession * gainSession = new G4UIGainServer();
  // "gainSession" is started.
  // Usage: gainSession->SessionStart();
  // "gainSession"  is deleted.
  // Usage: delete gainSession;
  //


   G4UIsession* SessionStart();  
   void PauseSessionStart(G4String msg);
   G4int ReceiveG4cout(G4String coutString);
   G4int ReceiveG4cerr(G4String cerrString);



public:
  //  G4UIGainServer(G4VUIshell* aShell=0);
   G4UIGainServer();
  //G4UIGainServer(void){}
    ~G4UIGainServer();
  //~G4UIGainServer(void){}

  //  void SetPrompt(const G4String& prompt);
  void SessionTerminate();
  void Prompt(G4String);
  G4String GetCommand();

private:
  void ExecuteCommand(G4String aCommand);
  G4bool GetHelpChoice(G4int& aInt);
  void ExitHelp();
  bool SetUPServer();
  void WaitingConnection();
  void CloseConnection();

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

  
};

#endif


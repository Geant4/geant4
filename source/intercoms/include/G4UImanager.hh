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
// $Id: G4UImanager.hh 102561 2017-02-09 08:16:05Z gcosmo $
//

#ifndef G4UImanager_h
#define G4UImanager_h 1

#include "globals.hh"

#include <vector>
#include <fstream>
#include "G4VStateDependent.hh"
#include "G4UIcommandStatus.hh"
class G4UIcommandTree;
class G4UIcommand;
class G4UIsession;
class G4UIcontrolMessenger;
class G4UnitsMessenger;
class G4LocalThreadCoutMessenger;
class G4UIaliasList;
class G4MTcoutDestination;
class G4UIbridge;

// class description:
//
//  This is a singlton class which controls the command manipulation
// and the user interface(s). The constructor of this class MUST NOT
// invoked by the user.
//

class G4UImanager : public G4VStateDependent
{
  public: // with description
      static G4UImanager * GetUIpointer();
      static G4UImanager * GetMasterUIpointer();
      //  A static method to get the pointer to the only existing object
      // of this class.

  protected:
      G4UImanager();
  public:
      ~G4UImanager();
  private:
      G4UImanager(const G4UImanager &right);
      const G4UImanager & operator=(const G4UImanager &right);
      G4int operator==(const G4UImanager &right) const;
      G4int operator!=(const G4UImanager &right) const;

  public: // with description
      G4String GetCurrentValues(const char * aCommand);
      //  This method returns a string which represents the current value(s)
      // of the parameter(s) of the specified command. Null string will be
      // returned if the given command is not defined or the command does
      // not support the GetCurrentValues() method.
      void AddNewCommand(G4UIcommand * newCommand);
      //  This method register a new command.
      void RemoveCommand(G4UIcommand * aCommand);
      //  This command remove the registered command. After invokation of this
      // command, that particular command cannot be applied.
      void ExecuteMacroFile(const char * fileName);
      //  A macro file defined by the argument will be read by G4UIbatch object.
      void Loop(const char * macroFile,const char * variableName,
                   G4double initialValue,G4double finalValue,G4double stepSize=1.0);
      // Execute a macro file more than once with a loop counter.
      void Foreach(const char * macroFile,const char * variableName,
                                                        const char * candidates);
      // Execute a macro file more than once with an aliased variable which takes
      // a value in the candidate list.
      G4int ApplyCommand(const char * aCommand);
      G4int ApplyCommand(const G4String& aCommand);
      //  A command (and parameter(s)) given
      // by the method's argument will be applied. Zero will be returned in
      // case the command is successfully executed. Positive non-zero value
      // will be returned if the command couldn't be executed. The meaning of
      // this non-zero value is the following.
      //   The returned number : xyy
      //        x00 : G4CommandStatus.hh enumeration
      //         yy : the problematic parameter (first found)
      void StoreHistory(const char* fileName = "G4history.macro");
      void StoreHistory(G4bool historySwitch,
                        const char* fileName = "G4history.macro");
      //  The executed commands will be stored in the defined file. If
      // "historySwitch" is false, saving will be suspended.
      void ListCommands(const char* direc);
      //  All commands registored under the given directory will be listed to
      // G4cout.
      void SetAlias(const char * aliasLine);
      //  Define an alias. The first word of "aliasLine" string is the
      // alias name and the remaining word(s) is(are) string value
      // to be aliased.
      void RemoveAlias(const char * aliasName);
      //  Remove the defined alias.
      void ListAlias();
      //  Print all aliases.
      G4String SolveAlias(const char* aCmd);
      //  Convert a command string which contains alias(es).
      void CreateHTML(const char* dir = "/");
      //  Generate HTML files for defined UI commands

  private:
      void AddWorkerCommand(G4UIcommand * newCommand);
      void RemoveWorkerCommand(G4UIcommand * aCommand);

  public:
      void LoopS(const char* valueList);
      void ForeachS(const char* valueList);
      //  These methods are used by G4UIcontrolMessenger to use Loop() and Foreach() methods.
      virtual G4bool Notify(G4ApplicationState requestedState);
      //  This method is exclusively invoked by G4StateManager and the user
      // must not use this method.

  private:
      void PauseSession(const char* msg);
      void CreateMessenger();
      G4UIcommandTree* FindDirectory(const char* dirName);

  //public:
  // following three methods will be removed quite soon.
  //    void Interact();
  //    void Interact(const char * promptCharacters);

  private:
      static G4ThreadLocal G4UImanager * fUImanager;
      static G4ThreadLocal G4bool fUImanagerHasBeenKilled;
      static G4UImanager * fMasterUImanager;
      G4UIcommandTree * treeTop;
      G4UIsession * session;
      G4UIsession * g4UIWindow;
      G4UIcontrolMessenger * UImessenger;
      G4UnitsMessenger * UnitsMessenger;
      G4LocalThreadCoutMessenger * CoutMessenger;
      G4String savedParameters;
      G4UIcommand * savedCommand;
      G4int verboseLevel;
      std::ofstream historyFile;
      G4bool saveHistory;
      std::vector<G4String> histVec;
      G4UIaliasList* aliasList;
      G4int maxHistSize;
      G4bool pauseAtBeginOfEvent;
      G4bool pauseAtEndOfEvent;
      G4String searchPath;
      std::vector<G4String> searchDirs;

  public: // with description
      G4String GetCurrentStringValue(const char * aCommand,
      G4int parameterNumber=1, G4bool reGet=true);
      G4int GetCurrentIntValue(const char * aCommand,
      G4int parameterNumber=1, G4bool reGet=true);
      G4double GetCurrentDoubleValue(const char * aCommand,
      G4int parameterNumber=1, G4bool reGet=true);
      G4String GetCurrentStringValue(const char * aCommand,
      const char * aParameterName, G4bool reGet=true);
      G4int GetCurrentIntValue(const char * aCommand,
      const char * aParameterName, G4bool reGet=true);
      G4double GetCurrentDoubleValue(const char * aCommand,
      const char * aParameterName, G4bool reGet=true);
      //  These six methods returns the current value of a parameter of the
      // given command. For the first three methods, the ordering number of
      // the parameter (1 is the first parameter) can be given, whereas,
      // other three methods can give the parameter name.
      //  If "reGet" is true, actual request of returning the current value
      // will be sent to the corresponding messenger, while, if it is false,
      // the value stored in G4Umanager will be used. The later case is valid
      // for the sequential invokation for the same command.

      inline void SetPauseAtBeginOfEvent(G4bool vl)
      { pauseAtBeginOfEvent = vl; }
      inline G4bool GetPauseAtBeginOfEvent() const
      { return pauseAtBeginOfEvent; }
      inline void SetPauseAtEndOfEvent(G4bool vl)
      { pauseAtEndOfEvent = vl; }
      inline G4bool GetPauseAtEndOfEvent() const
      { return pauseAtEndOfEvent; }
      //  If the boolean flags are true, Pause() method of G4StateManager is invoked
      // at the very begining (before generating a G4Event object) or at the end of
      // each event. So that, in case a (G)UI session is defined, the user can interact.


  public:
      inline G4UIcommandTree * GetTree() const
      { return treeTop; }
      inline G4UIsession * GetSession() const
      { return session; }
      inline G4UIsession * GetG4UIWindow() const
      { return g4UIWindow; }
  public: // with description
      inline void SetSession(G4UIsession *const value)
      { session = value; }
      inline void SetG4UIWindow(G4UIsession *const value)
      { g4UIWindow = value; }
      //  This method defines the active (G)UI session.
     void SetCoutDestination(G4UIsession *const value);
     //  This method defines the destination of G4cout/G4cerr stream.
     // For usual cases, this method will be invoked by a concrete
     // (G)UI session class object and thus the user needs not to invoke
     // this method by him(her)self.

  public:
      inline void SetVerboseLevel(G4int val)
      { verboseLevel = val; }
      inline G4int GetVerboseLevel() const
      { return verboseLevel; }
      inline G4int GetNumberOfHistory() const
      { return histVec.size(); }
      inline G4String GetPreviousCommand(G4int i) const
      {
        G4String st;
        if(i>=0 && i<G4int(histVec.size()))
        { st = histVec[i]; }
        return st;
      }
      inline void SetMaxHistSize(G4int mx)
      { maxHistSize = mx; }
      inline G4int GetMaxHistSize() const
      { return maxHistSize; }

      inline void SetMacroSearchPath(const G4String& path)
      { searchPath = path; }
      inline const G4String& GetMacroSearchPath() const
      { return searchPath; }
      void ParseMacroSearchPath();
      G4String FindMacroPath(const G4String& fname) const;

  private:
      G4bool isMaster;
      std::vector<G4UIbridge*>* bridges; 
      G4bool ignoreCmdNotFound;
      G4bool stackCommandsForBroadcast;
      std::vector<G4String>* commandStack;

  public:
      inline void SetMasterUIManager(G4bool val)
      {
        isMaster = val;
        //ignoreCmdNotFound = val;
        stackCommandsForBroadcast = val;
        if(val&&!bridges)
        {
           bridges = new std::vector<G4UIbridge*>;
           fMasterUImanager = this;
        }
      }
      inline void SetIgnoreCmdNotFound(G4bool val)
      { ignoreCmdNotFound = val; }

      std::vector<G4String>* GetCommandStack();
      void RegisterBridge(G4UIbridge* brg);

      void SetUpForAThread(G4int tId);
      //Setups as before but for a non-worker thread (e.g. vis)
      void SetUpForSpecialThread(G4String aPrefix);
    
      inline G4int GetThreadID() const
      { return threadID; }

  private:
      G4int threadID;
      G4MTcoutDestination* threadCout;
      static G4int igThreadID;

  public:
      void SetCoutFileName(const G4String& fileN = "G4cout.txt", G4bool ifAppend = true);
      void SetCerrFileName(const G4String& fileN = "G4cerr.txt", G4bool ifAppend = true);
      void SetThreadPrefixString(const G4String& s = "W");
      void SetThreadUseBuffer(G4bool flg = true);
      void SetThreadIgnore(G4int tid = 0);
      void SetThreadIgnoreInit(G4bool flg = true);
      inline G4MTcoutDestination* GetThreadCout() {return threadCout;};
 
  private:
      static G4bool doublePrecisionStr;

  public:
      inline static void UseDoublePrecisionStr(G4bool val) { doublePrecisionStr = val; }
      inline static G4bool DoublePrecisionStr() { return doublePrecisionStr; }

};

#endif

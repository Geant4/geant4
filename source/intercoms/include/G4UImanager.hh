// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UImanager.hh,v 1.4 1999-12-15 14:50:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4UImanager_h
#define G4UImanager_h 1

#include "globals.hh"
#include "g4rw/ctoken.h"
#include "g4rw/tvordvec.h"
#include "g4std/fstream"
class G4UIcommandTree;
class G4UIcommand;
class G4UIsession;
class G4UIcontrolMessenger;
class G4UnitsMessenger;

// class description:
//
//  This is a singlton class which controls the command manipulation
// and the user interface(s). The constructor of this class MUST NOT 
// invoked by the user.
//

class G4UImanager 
{
  public: // with description
      static G4UImanager * GetUIpointer();
      //  A static method to get the pointer to the only existing object
      // of this class.

  protected:
      G4UImanager();
  public:
      ~G4UImanager();
  private:
      G4UImanager(const G4UImanager &right);
      const G4UImanager & operator=(const G4UImanager &right);
      int operator==(const G4UImanager &right) const;
      int operator!=(const G4UImanager &right) const;

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
      void ExecuteMacroFile(G4String fileName);
      //  A macro file defined by the argument will be read by G4UIbatch object.
      G4int ApplyCommand(char * aCommand);
      G4int ApplyCommand(G4String aCommand);
      //  These two methods are identical. A command (and parameter(s)) given
      // by the method's argument will be applied. Zero will be returned in 
      // case the command is successfully executed. Positive non-zero value
      // will be returned if the command couldn't be executed. The meaning of
      // this non-zero value is the following.
      //   The returned number : xyy
      //        x00 : G4CommandStatus.hh enumeration
      //         yy : the problematic parameter (first found)
      void StoreHistory(G4String fileName="G4history.macro");
      void StoreHistory(G4bool historySwitch,
             G4String fileName="G4history.macro");
      //  The executed commands will be stored in the defined file. If 
      // "historySwitch" is false, saving will be suspended.
      void PauseSession(G4String msg);
      //  This method is exclusively invoked by G4StateManager and the user
      // must not use this method.
      void ListCommands(G4String direc);
      //  All commands registored under the given directory will be listed to
      // G4cout.
 
  private:
      void CreateMessenger();
      G4UIcommandTree* FindDirectory(const char* dirName);

  public:
  // following three methods will be removed quite soon.
      void Interact();
      void Interact(char * promptCharacters);
      void Interact(G4String promptCharacters);

  private:
      static G4UImanager * fUImanager;
      static G4bool fUImanagerHasBeenKilled;
      G4UIcommandTree * treeTop;
      G4UIsession * session;
      G4UIcontrolMessenger * UImessenger;
      G4UnitsMessenger * UnitsMessenger;
      G4String savedParameters;
      G4UIcommand * savedCommand;
      G4int verboseLevel;
      G4std::ofstream historyFile;
      G4bool saveHistory;
      G4RWTValOrderedVector<G4String> histVec;

  public: // with description
      G4String GetCurrentStringValue(const char * aCommand, 
	    int parameterNumber=1, G4bool reGet=true);
      G4int GetCurrentIntValue(const char * aCommand, 
	    int parameterNumber=1, G4bool reGet=true);
      G4double GetCurrentDoubleValue(const char * aCommand,
	    int parameterNumber=1, G4bool reGet=true);
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

  public:
      inline G4UIcommandTree * GetTree() const
      { return treeTop; };
      inline G4UIsession * GetSession() const
      { return session; };
  public: // with description
      inline void SetSession(G4UIsession *const value)
      { session = value; };
      //  This method defines the active (G)UI session.
     void SetCoutDestination(G4UIsession *const value);
     //  This method defines the destination of G4cout/G4cerr stream.
     // For usual cases, this method will be invoked by a concrete
     // (G)UI session class object and thus the user needs not to invoke
     // this method by him(her)self.

  public:
      inline void SetVerboseLevel(G4int val)
      { verboseLevel = val; };
      inline G4int GetVerboseLevel() const
      { return verboseLevel; };
      inline G4int GetNumberOfHistory() const
      { return histVec.entries(); };
      inline G4String GetPreviousCommand(G4int i) const
      { 
        G4String st;
        if(i>=0 && i<histVec.entries())
        { st = histVec[i]; }
        return st;
      }

  // Old methods kept for backward compatibility
  //    inline G4UIcommandTree * GetTree() const
  //    { return treeTop; };
  //    inline G4UIsession * GetSession() const
  //    { return session; };
  //    inline void SetSession(G4UIsession *const  value)
  //    { session = value; };

};

#endif


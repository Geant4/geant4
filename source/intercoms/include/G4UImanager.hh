// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UImanager.hh,v 1.1 1999-01-07 16:09:23 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4UImanager_h
#define G4UImanager_h 1

#include "globals.hh"
#include <rw/ctoken.h>
#include <rw/tvordvec.h>
#include <fstream.h>
class G4UIcommandTree;
class G4UIcommand;
class G4UIsession;
class G4UIcontrolMessenger;
class G4UnitsMessenger;

class G4UImanager 
{
  public:
      static G4UImanager * GetUIpointer();

  protected:
      G4UImanager();
  public:
      ~G4UImanager();
  private:
      G4UImanager(const G4UImanager &right);
      const G4UImanager & operator=(const G4UImanager &right);
      int operator==(const G4UImanager &right) const;
      int operator!=(const G4UImanager &right) const;

  public:
      G4String GetCurrentValues(const char * aCommand);
      void AddNewCommand(G4UIcommand * newCommand);
      void RemoveCommand(G4UIcommand * aCommand);
      void ExecuteMacroFile(G4String fileName);
      G4int ApplyCommand(char * aCommand);
      G4int ApplyCommand(G4String aCommand);
      void StoreHistory(G4String fileName="G4history.macro");
      void StoreHistory(G4bool historySwitch,
             G4String fileName="G4history.macro");
      void PauseSession(G4String msg);
      void ListCommands(G4String direc);
 
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
      ofstream historyFile;
      G4bool saveHistory;
      RWTValOrderedVector<G4String> histVec;

  public:
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

  public:
      inline G4UIcommandTree * GetTree() const
      { return treeTop; };
      inline G4UIsession * GetSession() const
      { return session; };
      inline void SetSession(G4UIsession *const value)
      { session = value; };
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

  public:
     void SetCoutDestination(G4UIsession *const value);


  // Old methods kept for backward compatibility
  //    inline G4UIcommandTree * GetTree() const
  //    { return treeTop; };
  //    inline G4UIsession * GetSession() const
  //    { return session; };
  //    inline void SetSession(G4UIsession *const  value)
  //    { session = value; };

};

#endif


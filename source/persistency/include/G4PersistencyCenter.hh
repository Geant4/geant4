// $Id: G4PersistencyCenter.hh,v 1.1 2002-11-24 13:45:23 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4PersistencyCenter.hh
//
// History:
//   '01.08.10  Youhei Morita  Initial creation (with "fadsclass3")

#ifndef PERSISTENCY_CENTER_HH
#define PERSISTENCY_CENTER_HH 1

#include <string>
#include <vector>
#include <map>
#include "G4HCIOcatalog.hh"
#include "G4VHCIOentry.hh"
#include "G4DCIOcatalog.hh"
#include "G4VDCIOentry.hh"
#include "G4FileUtilities.hh"

// Forward Declaration to avoid circular dependencies.
class G4PersistencyManager;

typedef std::map<std::string, G4PersistencyManager*,std::less<std::string> > PMap;
typedef std::map<int, std::string, std::less<int> > ObjMap;
typedef std::map<std::string, std::string, std::less<std::string> > FileMap;

enum StoreMode { kOn, kOff, kRecycle };

typedef std::map<std::string, StoreMode, std::less<std::string> > StoreMap;
typedef std::map<std::string, bool, std::less<std::string> >      BoolMap;

// Forward Declarations:
class G4PersistencyCenterMessenger;

// Class Description:
//   Class to handle loading of the G4PersistencyManager.

class G4PersistencyCenter
{
    public: // With description
      G4PersistencyCenter();
      G4PersistencyCenter(const G4PersistencyCenter&);
      // Constructor

      ~G4PersistencyCenter();
      // Destructor

    public: // With description
      static G4PersistencyCenter* GetG4PersistencyCenter();
      // returns the pointer of singleton G4PersistencyCenter

      void SelectSystem(std::string systemName);
      // Select the persistency package

      const std::string CurrentSystem() { return f_currentSystemName; };
      // returns the current persistent package name

      void SetHepMCObjyReaderFile(std::string file);
      // Sets the name of HepMCObjyReader file name.  To be called by generator.

      std::string CurrentHepMCObjyReaderFile();
      // Sets the name of HepMCObjyReader file name.  To be called by generator.

      void SetStoreMode(std::string objName, StoreMode mode);
      // Sets the object store mode.  Modes are kOn, kOff or kRecycle.

      void SetRetrieveMode(std::string objName, bool mode);
      // Sets the object retrieve mode.  Modes are true or false.

      StoreMode CurrentStoreMode(std::string objName);
      // returns the current object store mode.

      bool CurrentRetrieveMode(std::string objName);
      // returns the current object store mode.

      bool SetWriteFile(std::string objName, std::string writeFileName);
      // Sets the output filename.

      bool SetReadFile(std::string objName, std::string readFileName);
      // Sets the input filename.

      std::string CurrentWriteFile(std::string objName);
      // returns the current output filename.

      std::string CurrentReadFile(std::string objName);
      // returns the current input filename.

      std::string CurrentObject(std::string file);
      // returns the current object type

      void AddHCIOmanager(std::string detName, std::string colName);
      // add a hits colleciton I/O manager to the catalog

      std::string CurrentHCIOmanager();
      // Returns a list of registered hits colleciton I/O managers

      void AddDCIOmanager(std::string detName);
      // add a digits colleciton I/O manager to the catalog

      std::string CurrentDCIOmanager();
      // Returns a list of registered digits colleciton I/O managers

      void PrintAll();
      // prints the current G4PersistencyCenter settings.

      G4PersistencyManager* CurrentG4PersistencyManager() { return f_currentManager; };
      // returns the pointer of the currnet G4PersistencyManager.

      void SetG4PersistencyManager(G4PersistencyManager* pm, std::string name);
      // returns the pointer of the currnet G4PersistencyManager.

      G4PersistencyManager* GetG4PersistencyManager(std::string nam);
      // returns the pointer of the currnet G4PersistencyManager with name.

      void RegisterG4PersistencyManager(G4PersistencyManager* pm);
      // registers the persistency manager to the runtime catalog.

      void DeleteG4PersistencyManager();
      // deletes the current G4PersistencyManager.

      void SetVerboseLevel(int v);
      // Set verbose level.

      int VerboseLevel() { return m_verbose; };
      // Return verbose level.

    private:
      std::string PadString(std::string name, unsigned int width);
      // truncate or pad a string up to the width.

    private:
      G4PersistencyCenterMessenger* f_G4PersistencyCenterMessenger;

    private:
      G4PersistencyCenterMessenger* f_theMessenger;
      static G4PersistencyCenter*   f_thePointer;
      G4PersistencyManager*         f_currentManager;
      std::string                 f_currentSystemName;
      PMap                        f_theCatalog;
      ObjMap                      f_wrObj;
      ObjMap                      f_rdObj;
      FileMap                     f_writeFileName;
      FileMap                     f_readFileName;
      StoreMap                    f_writeFileMode;
      BoolMap                     f_readFileMode;
      int                         m_verbose;
      G4FileUtilities               f_ut;

}; // End of class G4PersistencyCenter

#endif


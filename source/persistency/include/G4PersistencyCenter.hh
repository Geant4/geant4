// $Id: G4PersistencyCenter.hh,v 1.3 2002-12-04 10:37:20 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4PersistencyCenter.hh
//
// History:
//   '01.08.10  Youhei Morita  Initial creation (with "fadsclass3")

#ifndef PERSISTENCY_CENTER_HH
#define PERSISTENCY_CENTER_HH 1

#include "G4Types.hh"
#include <string>
#include "g4std/vector"
#include "g4std/map"
#include "G4HCIOcatalog.hh"
#include "G4VHCIOentry.hh"
#include "G4DCIOcatalog.hh"
#include "G4VDCIOentry.hh"

#ifndef WIN32
  #include "G4FileUtilities.hh"
#endif

// Forward Declaration to avoid circular dependencies.
class G4PersistencyManager;

typedef G4std::map<G4std::string, G4PersistencyManager*,G4std::less<G4std::string> > PMap;
typedef G4std::map<int, G4std::string, G4std::less<int> > ObjMap;
typedef G4std::map<G4std::string, G4std::string, G4std::less<G4std::string> > FileMap;

enum StoreMode { kOn, kOff, kRecycle };

typedef G4std::map<G4std::string, StoreMode, G4std::less<G4std::string> > StoreMap;
typedef G4std::map<G4std::string, G4bool, G4std::less<G4std::string> >     BoolMap;

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

      void SelectSystem(G4std::string systemName);
      // Select the persistency package

      const G4std::string CurrentSystem() { return f_currentSystemName; };
      // returns the current persistent package name

      void SetHepMCObjyReaderFile(G4std::string file);
      // Sets the name of HepMCObjyReader file name.  To be called by generator.

      G4std::string CurrentHepMCObjyReaderFile();
      // Sets the name of HepMCObjyReader file name.  To be called by generator.

      void SetStoreMode(G4std::string objName, StoreMode mode);
      // Sets the object store mode.  Modes are kOn, kOff or kRecycle.

      void SetRetrieveMode(G4std::string objName, G4bool mode);
      // Sets the object retrieve mode.  Modes are true or false.

      StoreMode CurrentStoreMode(G4std::string objName);
      // returns the current object store mode.

      G4bool CurrentRetrieveMode(G4std::string objName);
      // returns the current object store mode.

      G4bool SetWriteFile(G4std::string objName, G4std::string writeFileName);
      // Sets the output filename.

      G4bool SetReadFile(G4std::string objName, G4std::string readFileName);
      // Sets the input filename.

      G4std::string CurrentWriteFile(G4std::string objName);
      // returns the current output filename.

      G4std::string CurrentReadFile(G4std::string objName);
      // returns the current input filename.

      G4std::string CurrentObject(G4std::string file);
      // returns the current object type

      void AddHCIOmanager(G4std::string detName, G4std::string colName);
      // add a hits colleciton I/O manager to the catalog

      G4std::string CurrentHCIOmanager();
      // Returns a list of registered hits colleciton I/O managers

      void AddDCIOmanager(G4std::string detName);
      // add a digits colleciton I/O manager to the catalog

      G4std::string CurrentDCIOmanager();
      // Returns a list of registered digits colleciton I/O managers

      void PrintAll();
      // prints the current G4PersistencyCenter settings.

      G4PersistencyManager* CurrentG4PersistencyManager() { return f_currentManager; };
      // returns the pointer of the currnet G4PersistencyManager.

      void SetG4PersistencyManager(G4PersistencyManager* pm, G4std::string name);
      // returns the pointer of the currnet G4PersistencyManager.

      G4PersistencyManager* GetG4PersistencyManager(G4std::string nam);
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
      G4std::string PadString(G4std::string name, unsigned int width);
      // truncate or pad a string up to the width.

    private:
      G4PersistencyCenterMessenger* f_G4PersistencyCenterMessenger;

    private:
      G4PersistencyCenterMessenger* f_theMessenger;
      static G4PersistencyCenter*   f_thePointer;
      G4PersistencyManager*         f_currentManager;
      G4std::string                 f_currentSystemName;
      PMap                          f_theCatalog;
      ObjMap                        f_wrObj;
      ObjMap                        f_rdObj;
      FileMap                       f_writeFileName;
      FileMap                       f_readFileName;
      StoreMap                      f_writeFileMode;
      BoolMap                       f_readFileMode;
      G4int                         m_verbose;
#ifndef WIN32
      G4FileUtilities               f_ut;
#endif
}; // End of class G4PersistencyCenter

#endif


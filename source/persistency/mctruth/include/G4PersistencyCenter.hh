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
// File: G4PersistencyCenter.hh
//
// History:
//   '01.08.10  Youhei Morita  Initial creation (with "fadsclass3")

#ifndef PERSISTENCY_CENTER_HH
#define PERSISTENCY_CENTER_HH 1

#include "G4Types.hh"
#include <string>
#include <vector>
#include <map>
#include "G4HCIOcatalog.hh"
#include "G4DCIOcatalog.hh"

#ifndef WIN32
  #include "G4FileUtilities.hh"
#endif

// Forward Declaration to avoid circular dependencies.
class G4PersistencyManager;

typedef std::map<std::string, G4PersistencyManager*,std::less<std::string> > PMap;
typedef std::map<int, std::string, std::less<int> > ObjMap;
typedef std::map<std::string, std::string, std::less<std::string> > FileMap;

enum StoreMode { kOn, kOff, kRecycle };

typedef std::map<std::string, StoreMode, std::less<std::string> > StoreMap;
typedef std::map<std::string, G4bool, std::less<std::string> >     BoolMap;

// Forward Declarations:
class G4PersistencyCenterMessenger;

// Class Description:
//   Class to handle loading of the G4PersistencyManager.

class G4PersistencyCenter
{
    public: // With description

      static G4PersistencyCenter* GetPersistencyCenter();
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

      void SetRetrieveMode(std::string objName, G4bool mode);
      // Sets the object retrieve mode.  Modes are true or false.

      StoreMode CurrentStoreMode(std::string objName);
      // returns the current object store mode.

      G4bool CurrentRetrieveMode(std::string objName);
      // returns the current object store mode.

      G4bool SetWriteFile(std::string objName, std::string writeFileName);
      // Sets the output filename.

      G4bool SetReadFile(std::string objName, std::string readFileName);
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

      G4PersistencyManager* CurrentPersistencyManager() { return f_currentManager; };
      // returns the pointer of the currnet G4PersistencyManager.

      void SetPersistencyManager(G4PersistencyManager* pm, std::string name);
      // returns the pointer of the currnet G4PersistencyManager.

      G4PersistencyManager* GetPersistencyManager(std::string nam);
      // returns the pointer of the currnet G4PersistencyManager with name.

      void RegisterPersistencyManager(G4PersistencyManager* pm);
      // registers the persistency manager to the runtime catalog.

      void DeletePersistencyManager();
      // deletes the current G4PersistencyManager.

      void SetVerboseLevel(int v);
      // Set verbose level.

      int VerboseLevel() { return m_verbose; };
      // Return verbose level.

    private:

      G4PersistencyCenter();
      // Constructor

      ~G4PersistencyCenter();
      // Destructor

      std::string PadString(std::string name, unsigned int width);
      // truncate or pad a string up to the width.

    private:
      G4PersistencyCenterMessenger* f_theMessenger;
      static G4ThreadLocal G4PersistencyCenter* f_thePointer;
      G4PersistencyManager*         f_currentManager;
      std::string                 f_currentSystemName;
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


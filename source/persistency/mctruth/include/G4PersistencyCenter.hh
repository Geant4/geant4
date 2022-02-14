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
// G4PersistencyCenter
//
// Class Description:
//
// Class to handle loading of the G4PersistencyManager.

// Author: Youhei Morita, 10.08.2001
// --------------------------------------------------------------------
#ifndef G4PERSISTENCYCENTER_HH
#define G4PERSISTENCYCENTER_HH 1

#include "G4Types.hh"
#include <string>
#include <vector>
#include <map>
#include "G4HCIOcatalog.hh"
#include "G4DCIOcatalog.hh"

#ifndef WIN32
#  include "G4FileUtilities.hh"
#endif

// Forward Declarations
class G4PersistencyManager;
class G4PersistencyCenterMessenger;

using PMap = std::map<G4String, G4PersistencyManager*, std::less<G4String>>;
using ObjMap = std::map<G4int, G4String, std::less<G4int>>;
using FileMap = std::map<G4String, G4String, std::less<G4String>>;

enum StoreMode
{
  kOn,
  kOff,
  kRecycle
};

using StoreMap = std::map<G4String, StoreMode, std::less<G4String>>;
using BoolMap = std::map<G4String, G4bool, std::less<G4String>>;

class G4PersistencyCenter
{
  public:

    static G4PersistencyCenter* GetPersistencyCenter();
      // Returns the pointer of singleton G4PersistencyCenter

    void SelectSystem(const G4String& systemName);
      // Select the persistency package

    const G4String& CurrentSystem() { return f_currentSystemName; }
      // Returns the current persistent package name

    void SetHepMCObjyReaderFile(const G4String& file);
      // Sets the name of HepMCObjyReader file name. To be called by generator

    G4String CurrentHepMCObjyReaderFile();
      // Sets the name of HepMCObjyReader file name. To be called by generator

    void SetStoreMode(const G4String& objName, StoreMode mode);
      // Sets the object store mode.  Modes are kOn, kOff or kRecycle

    void SetRetrieveMode(const G4String& objName, G4bool mode);
      // Sets the object retrieve mode.  Modes are true or false

    StoreMode CurrentStoreMode(const G4String& objName);
      // Returns the current object store mode

    G4bool CurrentRetrieveMode(const G4String& objName);
      // Returns the current object store mode

    G4bool SetWriteFile(const G4String& objName, const G4String& writeFileName);
      // Sets the output filename

    G4bool SetReadFile(const G4String& objName, const G4String& readFileName);
      // Sets the input filename

    G4String CurrentWriteFile(const G4String& objName);
      // Returns the current output filename

    G4String CurrentReadFile(const G4String& objName);
      // Returns the current input filename

    G4String CurrentObject(const G4String& file);
      // Returns the current object type

    void AddHCIOmanager(const G4String& detName, const G4String& colName);
      // Adds a hits colleciton I/O manager to the catalog

    G4String CurrentHCIOmanager();
      // Returns a list of registered hits colleciton I/O managers

    void AddDCIOmanager(const G4String& detName);
      // Adds a digits colleciton I/O manager to the catalog

    G4String CurrentDCIOmanager();
      // Returns a list of registered digits collection I/O managers

    void PrintAll();
      // Prints the current G4PersistencyCenter settings

    G4PersistencyManager* CurrentPersistencyManager()
    {
      return f_currentManager;
    }
      // Returns the pointer of the current G4PersistencyManager

    void SetPersistencyManager(G4PersistencyManager* pm, const G4String& name);
      // Returns the pointer of the current G4PersistencyManager

    G4PersistencyManager* GetPersistencyManager(const G4String& nam);
      // Returns the pointer of the current G4PersistencyManager with name

    void RegisterPersistencyManager(G4PersistencyManager* pm);
      // Registers the persistency manager to the runtime catalog

    void DeletePersistencyManager();
      // Deletes the current G4PersistencyManager

    void SetVerboseLevel(G4int v);
      // Sets verbose level

    G4int VerboseLevel() { return m_verbose; }
      // Returns verbose level

  private:

    G4PersistencyCenter();
      // Constructor

    ~G4PersistencyCenter();
      // Destructor

    G4String PadString(const G4String& name, unsigned int width);
      // Truncates or pad a string up to the width

  private:

    G4PersistencyCenterMessenger* f_theMessenger = nullptr;
    static G4ThreadLocal G4PersistencyCenter* f_thePointer;
    G4PersistencyManager* f_currentManager = nullptr;
    G4String f_currentSystemName;
    PMap f_theCatalog;
    ObjMap f_wrObj;
    ObjMap f_rdObj;
    FileMap f_writeFileName;
    FileMap f_readFileName;
    StoreMap f_writeFileMode;
    BoolMap f_readFileMode;
    G4int m_verbose = 0;
#ifndef WIN32
    G4FileUtilities f_ut;
#endif
};

#endif

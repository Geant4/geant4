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
// File: G4PersistencyCenter.cc
//
// History:
//   '01.08.10  Youhei Morita  Initial creation (with "fadsclass3")

#include "G4PersistencyCenter.hh"

// Forward Declarations:
#include "G4PersistencyCenterMessenger.hh"

// Addtional Include:
#include "G4UImanager.hh"
#include "G4PersistencyManager.hh"

G4PersistencyCenter* G4PersistencyCenter::f_thePointer=G4PersistencyCenter::GetPersistencyCenter();

// Implementation of Constructor #1
G4PersistencyCenter::G4PersistencyCenter()
 : m_verbose(0)
{
  f_wrObj[0] = "HepMC";
  f_wrObj[1] = "MCTruth";
  f_wrObj[2] = "Hits";
  f_wrObj[3] = "Digits";

  f_rdObj[0] = "Hits";
  f_rdObj[1] = "HitsBG";

  ObjMap::iterator itr;

  for ( itr = f_wrObj.begin(); itr != f_wrObj.end(); itr++ ) {
    f_writeFileName[(*itr).second] = "G4defaultOutput";
  }

  for ( itr = f_rdObj.begin(); itr != f_rdObj.end(); itr++ ) {
    f_readFileName[(*itr).second] = "G4defaultInput";
  }

  f_writeFileMode["HepMC"]   = kRecycle;
  f_writeFileMode["MCTruth"] = kOn;
  f_writeFileMode["Hits"]    = kOn;
  f_writeFileMode["Digits"]  = kOff;

  f_readFileMode["Hits"]   = false;
  f_readFileMode["HitsBG"] = false;

  f_theMessenger = new G4PersistencyCenterMessenger(this);
  f_currentManager = new G4PersistencyManager(this, "Default");
}

// Implementation of Constructor #2
G4PersistencyCenter::G4PersistencyCenter(const G4PersistencyCenter&)
{}

// Implementation of Destructor #1
G4PersistencyCenter::~G4PersistencyCenter()
{
  delete f_theMessenger;
  delete f_currentManager;
}

// Implementation of GetPersistencyCenter
G4PersistencyCenter* G4PersistencyCenter::GetPersistencyCenter()
{
  if ( f_thePointer == 0 ) f_thePointer = new G4PersistencyCenter;
  return f_thePointer;
}


// Implementation of SelectSystem
void G4PersistencyCenter::SelectSystem(G4std::string systemName)
{
  int st = 0;

  if (f_currentManager!=0) delete f_currentManager;

  G4PersistencyManager* pm = 0;

  if (systemName=="None")
  {
    G4cout<<" G4PersistencyCenter: Default is selected."<< G4endl;
    pm = new G4PersistencyManager(this, "Default");
  }
  else if (systemName=="ROOT")
  {
    G4cout<<" G4PersistencyCenter: \"ROOT\" Persistency Package is selected."
          <<G4endl;
    // G4UImanager *man=G4UImanager::GetUIpointer();
    // G4std::string libs="Cint:Core:Tree:Rint:Matrix:Physics:fadsROOT";
    // st = man->ApplyCommand("/load "+libs);
    if ( st == 0 ) {
      pm = GetPersistencyManager("ROOT");
    }
  }
  else if (systemName=="ODBMS")
  {
    G4cout<<" G4PersistencyCenter: \"ODBMS\" package is selected."<<G4endl;
    // G4UImanager *man=G4UImanager::GetUIpointer();
    // G4std::string libs="fadsODBMS";
    // st = man->ApplyCommand("/load "+libs);
    if ( st == 0 ) {
      pm = GetPersistencyManager("ODBMS");
    }
  }
  if ( st == 0 ) {
    f_currentManager = pm->Create();
    if (f_currentManager!=0) f_currentManager->SetVerboseLevel(m_verbose);
    f_currentSystemName = systemName;
  }
}

// Implementation of SetHepMCObjyReaderFile
void G4PersistencyCenter::SetHepMCObjyReaderFile(G4std::string file)
{
  if ( SetReadFile("HepMC", file) ) {
    SetRetrieveMode("HepMC", true);
  }
}

// Implementation of CurrentHepMCObjyReaderFile
G4std::string G4PersistencyCenter::CurrentHepMCObjyReaderFile()
{
  if ( CurrentRetrieveMode("HepMC") ) {
    return CurrentReadFile("HepMC");
  } else {
    return "";
  }
}

// Implementation of SetStoreMode
void G4PersistencyCenter::SetStoreMode(G4std::string objName, StoreMode mode)
{
  if ( (*(f_writeFileName.find(objName))).second != "" ) {
    f_writeFileMode[objName] = mode;
  } else {
    G4cerr << "!! unknown object type " << objName << " for output."
           << G4endl;
  }
}

// Implementation of SetRetrieveMode
void G4PersistencyCenter::SetRetrieveMode(G4std::string objName, G4bool mode)
{
  if ( (*(f_readFileName.find(objName))).second != "" ) {
    f_readFileMode[objName] = mode;
  } else {
    G4cerr << "!! unknown object type " << objName << " for input."
           << G4endl;
  }
}

// Implementation of CurrentStoreMode
StoreMode G4PersistencyCenter::CurrentStoreMode(G4std::string objName)
{


  if ( (*(f_writeFileName.find(objName))).second != "" ) {
    return f_writeFileMode[objName];
  } else {
    return kOff;
  }
}

// Implementation of CurrentRetrieveMode
G4bool G4PersistencyCenter::CurrentRetrieveMode(G4std::string objName)
{
  if ( (*(f_readFileName.find(objName))).second != "" ) {
    return f_readFileMode[objName];
  } else {
    return false;
  }
}

// Implementation of SetWriteFile
G4bool G4PersistencyCenter::SetWriteFile(G4std::string objName, G4std::string writeFileName)
{
  if ( (*(f_writeFileName.find(objName))).second != "" ) {
    f_writeFileName[objName] = writeFileName;
  } else {
    G4cerr << "!! unknown object type " << objName << " for output."
           << G4endl;
    return false;
  }
  return true;
}

// Implementation of SetReadFile
G4bool G4PersistencyCenter::SetReadFile(G4std::string objName, G4std::string readFileName)
{
#ifndef WIN32
  if ( f_ut.FileExists(readFileName) )
  {
    f_readFileName[objName] = readFileName;
  }
  else
  {
    G4cerr << "!! File \"" << objName << "\" does not exist."
           << G4endl;
    return false;
  }
#endif
  return true;
}

// Implementation of CurrentWriteFile
G4std::string G4PersistencyCenter::CurrentWriteFile(G4std::string objName)
{
  if ( (*(f_writeFileName.find(objName))).second != "" ) {
    return f_writeFileName[objName];
  } else {
    return "?????";
  }
}

// Implementation of CurrentReadFile
G4std::string G4PersistencyCenter::CurrentReadFile(G4std::string objName)
{
  if ( (*(f_readFileName.find(objName))).second != "" ) {
    return f_readFileName[objName];
  } else {
    return "?????";
  }
}

// Implementation of CurrentObject
G4std::string G4PersistencyCenter::CurrentObject(G4std::string file)
{
  FileMap::iterator itr;
  for ( itr = f_readFileName.begin(); itr != f_readFileName.end(); itr++ ) {
    if ( file == (*itr).second ) return (*itr).first;
  }
  for ( itr = f_writeFileName.begin(); itr != f_writeFileName.end(); itr++ ) {
    if ( file == (*itr).second ) return (*itr).first;
  }
  return "?????";
}

// Implementation of AddHCIOmanager
void G4PersistencyCenter::AddHCIOmanager(G4std::string detName, G4std::string colName)
{
  G4HCIOcatalog* ioc = G4HCIOcatalog::GetHCIOcatalog();

  G4VHCIOentry* ioe = ioc->GetEntry(detName);
  if ( ioe != 0 ) {
    ioe->CreateHCIOmanager(detName, colName);
  } else {
    G4cerr << "Error! -- HCIO assignment failed for detector " << detName
              << ", collection " << colName << G4endl;
  }
}

// Implementation of CurrentHCIOmanager
G4std::string G4PersistencyCenter::CurrentHCIOmanager()
{
  G4HCIOcatalog* ioc = G4HCIOcatalog::GetHCIOcatalog();
  return ioc->CurrentHCIOmanager();
}

// Implementation of AddDCIOmanager
void G4PersistencyCenter::AddDCIOmanager(G4std::string detName)
{
  G4DCIOcatalog* ioc = G4DCIOcatalog::GetDCIOcatalog();

  G4std::string colName = "";
  G4VDCIOentry* ioe = ioc->GetEntry(detName);
  if ( ioe != 0 ) {
    ioe->CreateDCIOmanager(detName, colName);
  } else {
    G4cerr << "Error! -- DCIO assignment failed for detector " << detName
              << ", collection " << colName << G4endl;
  }
}

// Implementation of CurrentDCIOmanager
G4std::string G4PersistencyCenter::CurrentDCIOmanager()
{
  G4DCIOcatalog* ioc = G4DCIOcatalog::GetDCIOcatalog();
  return ioc->CurrentDCIOmanager();
}

// Implementation of PrintAll
void G4PersistencyCenter::PrintAll()
{
  G4cout << "Persistency Package: " << CurrentSystem() << G4endl;
  G4cout << G4endl;

  ObjMap::iterator itr;
  G4std::string name;
  G4std::string file;
  StoreMode   mode;

  G4cout << "Output object types and file names:" << G4endl;
  for ( itr = f_wrObj.begin(); itr != f_wrObj.end(); itr++ ) {
    name = (*itr).second;
    // disabled HepMC and MCTruth for now
    if ( name != "HepMC" && name != "MCTruth" ) {
      G4cout << "  Object: " << PadString(name, 9);
      mode = CurrentStoreMode(name);
      if ( mode == kOn ) {
        G4cout << " <on>    ";
      } else if ( mode == kOff ) {
        G4cout << " <off>   ";
      } else if ( mode == kRecycle ) {
        G4cout << "<recycle>";
      }
      file = CurrentWriteFile(name);
      if ( file == "" ) file = "   <N/A>";
      G4cout << " File: " << file << G4endl;
    }
  }
  G4cout << G4endl;

  G4cout << "Input object types and file names:" << G4endl;
  for ( itr = f_rdObj.begin(); itr != f_rdObj.end(); itr++ ) {
    name = (*itr).second;
    // disabled HepMC and MCTruth for now
    if ( name != "HepMC" && name != "MCTruth" ) {
      G4cout << "  Object: " << PadString(name, 9);
      if ( CurrentRetrieveMode(name) ) {
        G4cout << " <on>    ";
      } else {
        G4cout << " <off>   ";
      }
      file = CurrentReadFile(name);
      if ( file == "" ) file = "   <N/A>";
      G4cout << " File: " << CurrentReadFile(name) << G4endl;
    }
  }
  G4cout << G4endl;

  G4HCIOcatalog* hioc = G4HCIOcatalog::GetHCIOcatalog();
  if ( hioc != 0 ) {
    G4cout << "Hit IO Managers:" << G4endl;
    hioc->PrintEntries();
    hioc->PrintHCIOmanager();
    G4cout << G4endl;
  } else {
    G4cout << "Hit IO Manager catalog is not registered." << G4endl;
  }

  G4DCIOcatalog* dioc = G4DCIOcatalog::GetDCIOcatalog();
  if ( dioc != 0 ) {
    G4cout << "Digit IO Managers:" << G4endl;
    dioc->PrintEntries();
    dioc->PrintDCIOmanager();
    G4cout << G4endl;
  } else {
    G4cout << "Digit IO Manager catalog is not registered." << G4endl;
  }
}

// Implementation of SetPersistencyManager
void G4PersistencyCenter::SetPersistencyManager(G4PersistencyManager* pm,
                                                  G4std::string name)
{
  f_currentManager=pm;
  f_currentSystemName=name;
}

// Implementation of GetPersistencyManager
G4PersistencyManager* G4PersistencyCenter::GetPersistencyManager(G4std::string nam)
{
  if (f_theCatalog.find(nam)!=f_theCatalog.end())
    return f_theCatalog[nam];
  return 0;
}

// Implementation of RegisterPersistencyManager
void G4PersistencyCenter::RegisterPersistencyManager(G4PersistencyManager* pm)
{
  f_theCatalog[pm->GetName()]=pm;
}

// Implementation of DeletePersistencyManager
void G4PersistencyCenter::DeletePersistencyManager()
{
  if (f_currentManager!=0) delete f_currentManager;
  f_currentManager=0;
}

// Implementation of SetVerboseLevel
void G4PersistencyCenter::SetVerboseLevel(int v)
{
  m_verbose = v;
  if ( f_currentManager != 0 ) f_currentManager->SetVerboseLevel(m_verbose);
}

// Implementation of PadString
G4std::string G4PersistencyCenter::PadString(G4std::string name, unsigned int width)
{
  if ( name.length() > width ) {
    return name.substr(0,width-1) + "#";
  } else {
    G4std::string wname = name;
    for ( unsigned int i=0; i < width-name.length(); i++) wname = wname + " ";
    return wname;
  }
}

// End of G4PersistencyCenter.cc


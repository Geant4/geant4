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
// G4PersistencyCenter implementation
//
// Author: Youhei Morita, 10.08.2001
// --------------------------------------------------------------------

#include "G4PersistencyCenter.hh"
#include "G4PersistencyCenterMessenger.hh"

#include "G4UImanager.hh"
#include "G4PersistencyManager.hh"
#include "G4VHCIOentry.hh"
#include "G4VDCIOentry.hh"

G4ThreadLocal G4PersistencyCenter* G4PersistencyCenter::f_thePointer = nullptr;

// --------------------------------------------------------------------
G4PersistencyCenter::G4PersistencyCenter()
{
  f_wrObj[0] = "HepMC";
  f_wrObj[1] = "MCTruth";
  f_wrObj[2] = "Hits";
  f_wrObj[3] = "Digits";

  f_rdObj[0] = "Hits";
  f_rdObj[1] = "HitsBG";

  for(auto itr = f_wrObj.cbegin(); itr != f_wrObj.cend(); ++itr)
  {
    f_writeFileName[(*itr).second] = "G4defaultOutput";
  }

  for(auto itr = f_rdObj.cbegin(); itr != f_rdObj.cend(); ++itr)
  {
    f_readFileName[(*itr).second] = "G4defaultInput";
  }

  f_writeFileMode["HepMC"]   = kRecycle;
  f_writeFileMode["MCTruth"] = kOn;
  f_writeFileMode["Hits"]    = kOn;
  f_writeFileMode["Digits"]  = kOff;

  f_readFileMode["Hits"]   = false;
  f_readFileMode["HitsBG"] = false;

  f_theMessenger   = new G4PersistencyCenterMessenger(this);
  f_currentManager = new G4PersistencyManager(this, "Default");
}

// --------------------------------------------------------------------
G4PersistencyCenter::~G4PersistencyCenter()
{
  delete f_theMessenger;
  delete f_currentManager;
}

// --------------------------------------------------------------------
G4PersistencyCenter* G4PersistencyCenter::GetPersistencyCenter()
{
  if(f_thePointer == nullptr)
    f_thePointer = new G4PersistencyCenter;
  return f_thePointer;
}

// --------------------------------------------------------------------
void G4PersistencyCenter::SelectSystem(const G4String& systemName)
{
  G4int st = 0;

  if(f_currentManager != nullptr)
    delete f_currentManager;

  G4PersistencyManager* pm = nullptr;

  if(systemName == "ROOT")
  {
    G4cout << " G4PersistencyCenter: \"ROOT\" Persistency Package is selected."
           << G4endl;
    // G4UImanager *man=G4UImanager::GetUIpointer();
    // std::string libs="Cint:Core:Tree:Rint:Matrix:Physics:fadsROOT";
    // st = man->ApplyCommand("/load "+libs);
    if(st == 0)
    {
      pm = GetPersistencyManager("ROOT");
    }
  }
  else if(systemName == "ODBMS")
  {
    G4cout << " G4PersistencyCenter: \"ODBMS\" package is selected." << G4endl;
    // G4UImanager *man=G4UImanager::GetUIpointer();
    // std::string libs="fadsODBMS";
    // st = man->ApplyCommand("/load "+libs);
    if(st == 0)
    {
      pm = GetPersistencyManager("ODBMS");
    }
  }
  else
  {
    G4cout << " G4PersistencyCenter: Default is selected." << G4endl;
    pm = new G4PersistencyManager(this, "Default");
  }

  if(st == 0)
  {
    f_currentManager = pm->Create();
    if(f_currentManager != nullptr)
      f_currentManager->SetVerboseLevel(m_verbose);
    f_currentSystemName = systemName;
  }
}

// --------------------------------------------------------------------
void G4PersistencyCenter::SetHepMCObjyReaderFile(const G4String& file)
{
  if(SetReadFile("HepMC", file))
  {
    SetRetrieveMode("HepMC", true);
  }
}

// --------------------------------------------------------------------
G4String G4PersistencyCenter::CurrentHepMCObjyReaderFile()
{
  if(CurrentRetrieveMode("HepMC"))
  {
    return CurrentReadFile("HepMC");
  }
  else
  {
    return "";
  }
}

// --------------------------------------------------------------------
void G4PersistencyCenter::SetStoreMode(const G4String& objName, StoreMode mode)
{
  if((*(f_writeFileName.find(objName))).second != "")
  {
    f_writeFileMode[objName] = mode;
  }
  else
  {
    G4cerr << "!! unknown object type " << objName << " for output." << G4endl;
  }
}

// --------------------------------------------------------------------
void G4PersistencyCenter::SetRetrieveMode(const G4String& objName, G4bool mode)
{
  if((*(f_readFileName.find(objName))).second != "")
  {
    f_readFileMode[objName] = mode;
  }
  else
  {
    G4cerr << "!! unknown object type " << objName << " for input." << G4endl;
  }
}

// --------------------------------------------------------------------
StoreMode G4PersistencyCenter::CurrentStoreMode(const G4String& objName)
{
  if((*(f_writeFileName.find(objName))).second != "")
  {
    return f_writeFileMode[objName];
  }
  else
  {
    return kOff;
  }
}

// --------------------------------------------------------------------
G4bool G4PersistencyCenter::CurrentRetrieveMode(const G4String& objName)
{
  if((*(f_readFileName.find(objName))).second != "")
  {
    return f_readFileMode[objName];
  }
  else
  {
    return false;
  }
}

// --------------------------------------------------------------------
G4bool G4PersistencyCenter::SetWriteFile(const G4String& objName,
                                         const G4String& writeFileName)
{
  if((*(f_writeFileName.find(objName))).second != "")
  {
    f_writeFileName[objName] = writeFileName;
  }
  else
  {
    G4cerr << "!! unknown object type " << objName << " for output." << G4endl;
    return false;
  }
  return true;
}

// --------------------------------------------------------------------
G4bool G4PersistencyCenter::SetReadFile(const G4String& objName,
                                        const G4String& readFileName)
{
#ifndef WIN32
  if(f_ut.FileExists(readFileName))
  {
    f_readFileName[objName] = readFileName;
  }
  else
  {
    G4cerr << "!! File \"" << objName << "\" does not exist." << G4endl;
    return false;
  }
#endif
  return true;
}

// --------------------------------------------------------------------
G4String G4PersistencyCenter::CurrentWriteFile(const G4String& objName)
{
  if((*(f_writeFileName.find(objName))).second != "")
  {
    return f_writeFileName[objName];
  }
  else
  {
    return "?????";
  }
}

// --------------------------------------------------------------------
G4String G4PersistencyCenter::CurrentReadFile(const G4String& objName)
{
  if((*(f_readFileName.find(objName))).second != "")
  {
    return f_readFileName[objName];
  }
  else
  {
    return "?????";
  }
}

// --------------------------------------------------------------------
G4String G4PersistencyCenter::CurrentObject(const G4String& file)
{
  for(auto itr = f_readFileName.cbegin(); itr != f_readFileName.cend(); ++itr)
  {
    if(file == (*itr).second)
      return (*itr).first;
  }
  for(auto itr = f_writeFileName.cbegin(); itr != f_writeFileName.cend(); ++itr)
  {
    if(file == (*itr).second)
      return (*itr).first;
  }
  return "?????";
}

// --------------------------------------------------------------------
void G4PersistencyCenter::AddHCIOmanager(const G4String& detName,
                                         const G4String& colName)
{
  G4HCIOcatalog* ioc = G4HCIOcatalog::GetHCIOcatalog();

  G4VHCIOentry* ioe = ioc->GetEntry(detName);
  if(ioe != nullptr)
  {
    ioe->CreateHCIOmanager(detName, colName);
  }
  else
  {
    G4cerr << "Error! -- HCIO assignment failed for detector " << detName
           << ", collection " << colName << G4endl;
  }
}

// --------------------------------------------------------------------
G4String G4PersistencyCenter::CurrentHCIOmanager()
{
  G4HCIOcatalog* ioc = G4HCIOcatalog::GetHCIOcatalog();
  return ioc->CurrentHCIOmanager();
}

// --------------------------------------------------------------------
void G4PersistencyCenter::AddDCIOmanager(const G4String& detName)
{
  G4DCIOcatalog* ioc = G4DCIOcatalog::GetDCIOcatalog();

  G4String colName = "";
  G4VDCIOentry* ioe = ioc->GetEntry(detName);
  if(ioe != nullptr)
  {
    ioe->CreateDCIOmanager(detName, colName);
  }
  else
  {
    G4cerr << "Error! -- DCIO assignment failed for detector " << detName
           << ", collection " << colName << G4endl;
  }
}

// --------------------------------------------------------------------
G4String G4PersistencyCenter::CurrentDCIOmanager()
{
  G4DCIOcatalog* ioc = G4DCIOcatalog::GetDCIOcatalog();
  return ioc->CurrentDCIOmanager();
}

// --------------------------------------------------------------------
void G4PersistencyCenter::PrintAll()
{
  G4cout << "Persistency Package: " << CurrentSystem() << G4endl;
  G4cout << G4endl;

  G4String name, file;
  StoreMode mode;

  G4cout << "Output object types and file names:" << G4endl;
  for(auto itr = f_wrObj.cbegin(); itr != f_wrObj.cend(); ++itr)
  {
    name = (*itr).second;
    // disabled HepMC and MCTruth for now
    if(name != "HepMC" && name != "MCTruth")
    {
      G4cout << "  Object: " << PadString(name, 9);
      mode = CurrentStoreMode(name);
      if(mode == kOn)
      {
        G4cout << " <on>    ";
      }
      else if(mode == kOff)
      {
        G4cout << " <off>   ";
      }
      else if(mode == kRecycle)
      {
        G4cout << "<recycle>";
      }
      file = CurrentWriteFile(name);
      if(file == "")
        file = "   <N/A>";
      G4cout << " File: " << file << G4endl;
    }
  }
  G4cout << G4endl;

  G4cout << "Input object types and file names:" << G4endl;
  for(auto itr = f_rdObj.cbegin(); itr != f_rdObj.cend(); ++itr)
  {
    name = (*itr).second;
    // disabled HepMC and MCTruth for now
    if(name != "HepMC" && name != "MCTruth")
    {
      G4cout << "  Object: " << PadString(name, 9);
      if(CurrentRetrieveMode(name))
      {
        G4cout << " <on>    ";
      }
      else
      {
        G4cout << " <off>   ";
      }
      file = CurrentReadFile(name);
      if(file == "")
        file = "   <N/A>";
      G4cout << " File: " << CurrentReadFile(name) << G4endl;
    }
  }
  G4cout << G4endl;

  G4HCIOcatalog* hioc = G4HCIOcatalog::GetHCIOcatalog();
  if(hioc != nullptr)
  {
    G4cout << "Hit IO Managers:" << G4endl;
    hioc->PrintEntries();
    hioc->PrintHCIOmanager();
    G4cout << G4endl;
  }
  else
  {
    G4cout << "Hit IO Manager catalog is not registered." << G4endl;
  }

  G4DCIOcatalog* dioc = G4DCIOcatalog::GetDCIOcatalog();
  if(dioc != nullptr)
  {
    G4cout << "Digit IO Managers:" << G4endl;
    dioc->PrintEntries();
    dioc->PrintDCIOmanager();
    G4cout << G4endl;
  }
  else
  {
    G4cout << "Digit IO Manager catalog is not registered." << G4endl;
  }
}

// --------------------------------------------------------------------
void G4PersistencyCenter::SetPersistencyManager(G4PersistencyManager* pm,
                                                const G4String& name)
{
  f_currentManager    = pm;
  f_currentSystemName = name;
}

// --------------------------------------------------------------------
G4PersistencyManager*
G4PersistencyCenter::GetPersistencyManager(const G4String& nam)
{
  if(f_theCatalog.find(nam) != f_theCatalog.cend())
    return f_theCatalog[nam];
  return nullptr;
}

// --------------------------------------------------------------------
void G4PersistencyCenter::RegisterPersistencyManager(G4PersistencyManager* pm)
{
  f_theCatalog[pm->GetName()] = pm;
}

// --------------------------------------------------------------------
void G4PersistencyCenter::DeletePersistencyManager()
{
  if(f_currentManager != nullptr)
    delete f_currentManager;
  f_currentManager = nullptr;
}

// --------------------------------------------------------------------
void G4PersistencyCenter::SetVerboseLevel(G4int v)
{
  m_verbose = v;
  if(f_currentManager != nullptr)
    f_currentManager->SetVerboseLevel(m_verbose);
}

// --------------------------------------------------------------------
G4String G4PersistencyCenter::PadString(const G4String& name,
                                        unsigned int width)
{
  if(name.length() > width)
  {
    return name.substr(0, width - 1) + "#";
  }
  else
  {
    G4String wname = name;
    for(unsigned int i = 0; i < width - name.length(); ++i)
      wname = wname + " ";
    return wname;
  }
}

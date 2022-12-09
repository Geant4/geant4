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
// G4PersistencyCenterMessenger implementation
//
// Author: Youhei Morita, 18.07.2001
// --------------------------------------------------------------------

#include "G4PersistencyCenterMessenger.hh"

// --------------------------------------------------------------------
G4PersistencyCenterMessenger::
G4PersistencyCenterMessenger(G4PersistencyCenter* p)
  : pc(p)
{
  G4String name = "/persistency/";
  directory        = new G4UIdirectory(name.c_str());
  directory->SetGuidance("Control commands for Persistency package");

  G4String cmd = name + "verbose";

  verboseCmd = new G4UIcmdWithAnInteger(cmd.c_str(), this);
  verboseCmd->SetGuidance("Set the verbose level of G4PersistencyManager.");
  verboseCmd->SetGuidance(" 0 : Silent (default)");
  verboseCmd->SetGuidance(" 1 : Display main topics");
  verboseCmd->SetGuidance(" 2 : Display event-level topics");
  verboseCmd->SetGuidance(" 3 : Display debug information");
  verboseCmd->SetParameterName("level", true);
  verboseCmd->SetDefaultValue(0);
  verboseCmd->SetRange("level >=0 && level <=3");

  G4String vname = name + "select";

  cmd    = vname;
  select = new G4UIcmdWithAString(cmd.c_str(), this);
  select->SetGuidance("Selection of a persistency package");
  select->SetParameterName("Persistency package name", true, true);
  select->SetCandidates("ODBMS ROOT None");

  vname = name + "store/";

  subdir1 = new G4UIdirectory(vname.c_str());
  subdir1->SetGuidance("Specifiy object types for store");

  wrObj.push_back("HepMC");
  wrObj.push_back("MCTruth");
  wrObj.push_back("Hits");

  G4String guidance;
  G4int i;

  for(i = 0; i < 3; ++i)
  {
    cmd      = vname + wrObj[i];
    guidance = "Store " + wrObj[i] + " objects for output";
    storeObj.push_back(new G4UIcmdWithAString(cmd.c_str(), this));
    storeObj[i]->SetGuidance(guidance.c_str());
    if(wrObj[i] == "HepMC")
    {
      storeObj[i]->SetCandidates("on off recycle");
    }
    else
    {
      storeObj[i]->SetCandidates("on off");
    }
  }

  vname += "using/";
  subdir2 = new G4UIdirectory(vname.c_str());
  subdir2->SetGuidance("Select I/O manager for store");

  cmd      = vname + "hitIO";
  regHitIO = new G4UIcmdWithAString(cmd.c_str(), this);
  regHitIO->SetGuidance("Resiter Hits I/O Manager");
  regHitIO->SetParameterName("Name of Hits I/O Manager", true, true);

  vname   = name + "set/";
  subdir3 = new G4UIdirectory(vname.c_str());
  subdir3->SetGuidance("Set various parameters");

  vname += "writeFile/";
  subdir4 = new G4UIdirectory(vname.c_str());
  subdir4->SetGuidance("Set output file names for object types");

  for(i = 0; i < 3; ++i)
  {
    cmd      = vname + wrObj[i];
    guidance = "Set an output file name for " + wrObj[i] + ".";
    setWrFile.push_back(new G4UIcmdWithAString(cmd.c_str(), this));
    setWrFile[i]->SetGuidance(guidance.c_str());
    setWrFile[i]->SetParameterName("file name", true, true);
  }

  vname   = name + "set/ReadFile/";
  subdir5 = new G4UIdirectory(vname.c_str());
  subdir5->SetGuidance("Set input file names for object types");

  rdObj.push_back("Hits");

  cmd      = vname + rdObj[0];
  guidance = "Set an input file name for " + rdObj[0] + ".";
  setRdFile.push_back(new G4UIcmdWithAString(cmd.c_str(), this));
  setRdFile[0]->SetGuidance(guidance.c_str());
  setRdFile[0]->SetParameterName("file name", true, true);

  cmd      = name + "printall";
  printAll = new G4UIcmdWithoutParameter(cmd.c_str(), this);
  printAll->SetGuidance("Print all parameters.");
}

// --------------------------------------------------------------------
G4PersistencyCenterMessenger::~G4PersistencyCenterMessenger()
{
  delete directory;
  delete subdir1;
  delete subdir2;
  delete subdir3;
  delete subdir4;
  delete subdir5;
  delete verboseCmd;
  delete select;
  delete regHitIO;
  for(G4int i = 0; i < 3; ++i)
  {
    delete storeObj[i];
    delete setWrFile[i];
  }
  delete setRdFile[0];
  delete printAll;
}

// --------------------------------------------------------------------
void G4PersistencyCenterMessenger::SetNewValue(G4UIcommand* command,
                                               G4String newValues)
{
  if(command == verboseCmd)
  {
    pc->SetVerboseLevel(verboseCmd->GetNewIntValue(newValues));
  }
  else if(command == select)
  {
    pc->SelectSystem(newValues);
  }
  else if(command == regHitIO)
  {
    pc->AddHCIOmanager(PopWord(newValues, 1, " "), PopWord(newValues, 2, " "));
  }
  else if(command == setRdFile[0])
  {
    pc->SetReadFile(rdObj[0], newValues);
  }
  else if(command == printAll)
  {
    pc->PrintAll();
  }
  else
  {
    for(G4int i = 0; i < 3; ++i)
    {
      if(command == storeObj[i])
      {
        StoreMode mode = kOff;
        if(newValues == "on")
        {
          mode = kOn;
        }
        else if(newValues == "off")
        {
          mode = kOff;
        }
        else if(newValues == "recycle")
        {
          mode = kRecycle;
        }
        else
        {
          G4cerr << "Unrecognized keyword - \"" << newValues << "\"." << G4endl;
        }
        pc->SetStoreMode(wrObj[i], mode);
        break;
      }
      else if(command == setWrFile[i])
      {
        pc->SetWriteFile(wrObj[i], newValues);
        break;
      }
    }
  }
}

// --------------------------------------------------------------------
G4String G4PersistencyCenterMessenger::GetCurrentValue(G4UIcommand* command)
{
  G4String ustr = "Undefined";

  if(command == verboseCmd)
  {
    return G4UIcommand::ConvertToString(pc->VerboseLevel());
  }
  else if(command == select)
  {
    return pc->CurrentSystem();
  }
  else if(command == regHitIO)
  {
    return pc->CurrentHCIOmanager();
  }
  else if(command == setRdFile[0])
  {
    return pc->CurrentReadFile(rdObj[0]);
  }
  else
  {
    for(G4int i = 0; i < 3; ++i)
    {
      if(command == storeObj[i])
      {
        switch(pc->CurrentStoreMode(wrObj[i]))
        {
          case kOn:
            return "on";
            break;
          case kOff:
            return "off";
            break;
          case kRecycle:
            return "recycle";
            break;
          default:
            return "?????";
            break;
        };
      }
      else if(command == setWrFile[i])
      {
        return pc->CurrentWriteFile(wrObj[i]);
      }
    }
  }

  return ustr;
}

// --------------------------------------------------------------------
G4String G4PersistencyCenterMessenger::PopWord(const G4String& text, G4int n,
                                               const G4String& delim)
{
  if(text.length() <= 0)
    return "";
  std::size_t p = 0, p0 = 0;
  std::size_t p1 = 0;
  for(G4int i = 0; i < n; ++i)
  {
    p1 = text.find_first_of(delim, p0 + 1);
    while(p1 == p0 + 1)
    {
      p0 = p1;
      p1 = text.find_first_of(delim, p0 + 1);
    }
    p = p0;
    if(p1 == G4String::npos)
    {
      if(i + 1 < n)
        return "";
      p1 = text.length();
      break;
    }
    p0 = p1;
  }
  if(p > 0)
    ++p;
  return text.substr(p, p1 - p);
}

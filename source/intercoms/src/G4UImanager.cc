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
// G4UImanager
//
// Author: Makoto Asai, 1997
// --------------------------------------------------------------------

#include "G4UImanager.hh"
#include "G4UIcommandTree.hh"
#include "G4UIcommand.hh"
#include "G4UIsession.hh"
#include "G4UIbatch.hh"
#include "G4UIcontrolMessenger.hh"
#include "G4UnitsMessenger.hh"
#include "G4LocalThreadCoutMessenger.hh"
#include "G4ProfilerMessenger.hh"
#include "G4ios.hh"
#include "G4strstreambuf.hh"
#include "G4StateManager.hh"
#include "G4UIaliasList.hh"
#include "G4Tokenizer.hh"
#include "G4MTcoutDestination.hh"
#include "G4UIbridge.hh"
#include "G4Threading.hh"

#include <sstream>
#include <fstream>

G4bool G4UImanager::doublePrecisionStr = false;
G4int G4UImanager::igThreadID = -1;

// --------------------------------------------------------------------
G4UImanager*& G4UImanager::fUImanager()
{
  G4ThreadLocalStatic G4UImanager* _instance = nullptr;
  return _instance;
}

// --------------------------------------------------------------------
G4bool& G4UImanager::fUImanagerHasBeenKilled()
{
  G4ThreadLocalStatic G4bool _instance = false;
  return _instance;
}

// --------------------------------------------------------------------
G4UImanager*& G4UImanager::fMasterUImanager()
{
  static G4UImanager* _instance = nullptr;
  return _instance;
}

// --------------------------------------------------------------------
G4UImanager* G4UImanager::GetUIpointer()
{
  if(fUImanager() == nullptr)
  {
    if(!fUImanagerHasBeenKilled())
    {
      fUImanager() = new G4UImanager;
      fUImanager()->CreateMessenger();
    }
  }
  return fUImanager();
}

// --------------------------------------------------------------------
G4UImanager* G4UImanager::GetMasterUIpointer()
{
  return fMasterUImanager();
}

// --------------------------------------------------------------------
G4UImanager::G4UImanager()
  : G4VStateDependent(true)
{
  treeTop      = new G4UIcommandTree("/");
  aliasList    = new G4UIaliasList;
  SetCoutDestination(session);
  commandStack = new std::vector<G4String>;
}

// --------------------------------------------------------------------
void G4UImanager::CreateMessenger()
{
  UImessenger      = new G4UIcontrolMessenger;
  UnitsMessenger   = new G4UnitsMessenger;
  CoutMessenger    = new G4LocalThreadCoutMessenger;
  ProfileMessenger = new G4ProfilerMessenger;
}

// --------------------------------------------------------------------
G4UImanager::~G4UImanager()
{
  if(bridges != nullptr)
  {
    for(auto bridge : *bridges)
    {
      delete bridge;
    }
    delete bridges;
  }
  SetCoutDestination(nullptr);
  histVec.clear();
  if(saveHistory)
  {
    historyFile.close();
  }
  delete CoutMessenger;
  delete ProfileMessenger;
  delete UnitsMessenger;
  delete UImessenger;
  delete treeTop;
  delete aliasList;
  fUImanagerHasBeenKilled() = true;
  fUImanager()              = nullptr;
  if(commandStack != nullptr)
  {
    commandStack->clear();
    delete commandStack;
  }
  if(threadID >= 0)
  {
    delete threadCout;
    G4iosFinalization();
    threadID = -1;
  }
}

// --------------------------------------------------------------------
void G4UImanager::UseDoublePrecisionStr(G4bool val)
{
  doublePrecisionStr = val;
}

// --------------------------------------------------------------------
G4bool G4UImanager::DoublePrecisionStr()
{
  return doublePrecisionStr;
}

// --------------------------------------------------------------------
G4String G4UImanager::GetCurrentValues(const char* aCommand)
{
  G4String theCommand = aCommand;
  savedCommand        = treeTop->FindPath(theCommand);
  if(savedCommand == nullptr)
  {
    G4cerr << "command not found" << G4endl;
    return G4String();
  }
  return savedCommand->GetCurrentValue();
}

// --------------------------------------------------------------------
G4String G4UImanager::GetCurrentStringValue(const char* aCommand,
                                            G4int parameterNumber, G4bool reGet)
{
  if(reGet || savedCommand == nullptr)
  {
    savedParameters = GetCurrentValues(aCommand);
  }
  G4Tokenizer savedToken(savedParameters);
  G4String token;
  for(G4int i_thParameter = 0; i_thParameter < parameterNumber; ++i_thParameter)
  {
    token = savedToken();
    if(token.empty())
    {
      return G4String();
    }
    if(token[(size_t) 0] == '"')
    {
      token.append(" ");
      token.append(savedToken("\""));
    }
  }
  return token;
}

// --------------------------------------------------------------------
G4String G4UImanager::GetCurrentStringValue(const char* aCommand,
                                            const char* aParameterName,
                                            G4bool reGet)
{
  if(reGet || savedCommand == nullptr)
  {
    G4String parameterValues = GetCurrentValues(aCommand);
  }
  for(G4int i = 0; i < (G4int)savedCommand->GetParameterEntries(); ++i)
  {
    if(aParameterName == savedCommand->GetParameter(i)->GetParameterName())
    {
      return GetCurrentStringValue(aCommand, i + 1, false);
    }
  }
  return G4String();
}

// --------------------------------------------------------------------
G4int G4UImanager::GetCurrentIntValue(const char* aCommand,
                                      const char* aParameterName, G4bool reGet)
{
  G4String targetParameter =
    GetCurrentStringValue(aCommand, aParameterName, reGet);
  G4int value;
  const char* t = targetParameter;
  std::istringstream is(t);
  is >> value;
  return value;
}

// --------------------------------------------------------------------
G4int G4UImanager::GetCurrentIntValue(const char* aCommand,
                                      G4int parameterNumber, G4bool reGet)
{
  G4String targetParameter =
    GetCurrentStringValue(aCommand, parameterNumber, reGet);
  G4int value;
  const char* t = targetParameter;
  std::istringstream is(t);
  is >> value;
  return value;
}

// --------------------------------------------------------------------
G4double G4UImanager::GetCurrentDoubleValue(const char* aCommand,
                                            const char* aParameterName,
                                            G4bool reGet)
{
  G4String targetParameter =
    GetCurrentStringValue(aCommand, aParameterName, reGet);
  G4double value;
  const char* t = targetParameter;
  std::istringstream is(t);
  is >> value;
  return value;
}

// --------------------------------------------------------------------
G4double G4UImanager::GetCurrentDoubleValue(const char* aCommand,
                                            G4int parameterNumber, G4bool reGet)
{
  G4String targetParameter =
    GetCurrentStringValue(aCommand, parameterNumber, reGet);
  G4double value;
  const char* t = targetParameter;
  std::istringstream is(t);
  is >> value;
  return value;
}

// --------------------------------------------------------------------
void G4UImanager::AddNewCommand(G4UIcommand* newCommand)
{
  treeTop->AddNewCommand(newCommand);
  if(fMasterUImanager() != nullptr && G4Threading::G4GetThreadId() == 0)
  {
    fMasterUImanager()->AddWorkerCommand(newCommand);
  }
}

// --------------------------------------------------------------------
void G4UImanager::AddWorkerCommand(G4UIcommand* newCommand)
{
  treeTop->AddNewCommand(newCommand, true);
}

// --------------------------------------------------------------------
void G4UImanager::RemoveCommand(G4UIcommand* aCommand)
{
  treeTop->RemoveCommand(aCommand);
  if(fMasterUImanager() != nullptr && G4Threading::G4GetThreadId() == 0)
  {
    fMasterUImanager()->RemoveWorkerCommand(aCommand);
  }
}

// --------------------------------------------------------------------
void G4UImanager::RemoveWorkerCommand(G4UIcommand* aCommand)
{
  treeTop->RemoveCommand(aCommand, true);
}

// --------------------------------------------------------------------
void G4UImanager::ExecuteMacroFile(const char* fileName)
{
  G4UIsession* batchSession    = new G4UIbatch(fileName, session);
  session                      = batchSession;
  lastRC                       = 0;
  G4UIsession* previousSession = session->SessionStart();
  lastRC                       = session->GetLastReturnCode();
  delete session;
  session = previousSession;
}

// --------------------------------------------------------------------
void G4UImanager::LoopS(const char* valueList)
{
  G4String vl = valueList;
  G4Tokenizer parameterToken(vl);
  G4String mf = parameterToken();
  G4String vn = parameterToken();
  G4String c1 = parameterToken();
  c1 += " ";
  c1 += parameterToken();
  c1 += " ";
  c1 += parameterToken();
  const char* t1 = c1;
  std::istringstream is(t1);
  G4double d1;
  G4double d2;
  G4double d3;
  is >> d1 >> d2 >> d3;
  Loop(mf, vn, d1, d2, d3);
}

// --------------------------------------------------------------------
void G4UImanager::Loop(const char* macroFile, const char* variableName,
                       G4double initialValue, G4double finalValue,
                       G4double stepSize)
{
  G4String cd;
  if(stepSize > 0)
  {
    for(G4double d = initialValue; d <= finalValue; d += stepSize)
    {
      std::ostringstream os;
      os << d;
      cd += os.str();
      cd += " ";
    }
  }
  else
  {
    for(G4double d = initialValue; d >= finalValue; d += stepSize)
    {
      std::ostringstream os;
      os << d;
      cd += os.str();
      cd += " ";
    }
  }
  Foreach(macroFile, variableName, cd);
}

// --------------------------------------------------------------------
void G4UImanager::ForeachS(const char* valueList)
{
  G4String vl = valueList;
  G4Tokenizer parameterToken(vl);
  G4String mf = parameterToken();
  G4String vn = parameterToken();
  G4String c1 = parameterToken();
  G4String ca;
  while(!((ca = parameterToken()).empty()))
  {
    c1 += " ";
    c1 += ca;
  }

  G4String aliasValue = c1;
  if(aliasValue[0] == '"')
  {
    G4String strippedValue;
    if(aliasValue.back() == '"')
    {
      strippedValue = aliasValue.substr(1, aliasValue.length() - 2);
    }
    else
    {
      strippedValue = aliasValue.substr(1, aliasValue.length() - 1);
    }
    aliasValue = strippedValue;
  }

  //  Foreach(mf,vn,c1);
  Foreach(mf, vn, aliasValue);
}

// --------------------------------------------------------------------
void G4UImanager::Foreach(const char* macroFile, const char* variableName,
                          const char* candidates)
{
  G4String candidatesString = candidates;
  G4Tokenizer parameterToken(candidatesString);
  G4String cd;
  while(!((cd = parameterToken()).empty()))
  {
    G4String vl = variableName;
    vl += " ";
    vl += cd;
    SetAlias(vl);
    ExecuteMacroFile(FindMacroPath(macroFile));
    if(lastRC != 0)
    {
      G4ExceptionDescription ed;
      ed << "Loop aborted due to a command execution error - "
         << "error code " << lastRC;
      G4Exception("G4UImanager::Foreach", "UIMAN0201", JustWarning, ed);
      break;
    }
  }
}

// --------------------------------------------------------------------
G4String G4UImanager::SolveAlias(const char* aCmd)
{
  G4String aCommand = aCmd;
  std::size_t ia = aCommand.find("{");
  std::size_t iz = aCommand.find("#");
  while((ia != std::string::npos) &&
        ((iz == std::string::npos) || (ia < iz)))
  {
    G4int ibx = -1;
    while(ibx < 0)
    {
      std::size_t ib = aCommand.find("}");
      if(ib == std::string::npos)
      {
        G4cerr << aCommand << G4endl;
        for(std::size_t i = 0; i < ia; ++i)
        {
          G4cerr << " ";
        }
        G4cerr << "^" << G4endl;
        G4cerr << "Unmatched alias parenthesis -- command ignored" << G4endl;
        G4String nullStr;
        return nullStr;
      }
      G4String ps = aCommand.substr(ia + 1, aCommand.length() - (ia + 1));
      std::size_t ic    = ps.find("{");
      std::size_t id    = ps.find("}");
      if(ic != std::string::npos && ic < id)
      {
        ia += ic + 1;
      }
      else
      {
        ibx = (G4int)ib;
      }
    }
    //--- Here ia represents the position of innermost "{"
    //--- and ibx represents corresponding "}"
    G4String subs;
    if(ia > 0)
    {
      subs = aCommand.substr(0, ia);
    }
    G4String alis = aCommand.substr(ia + 1, ibx - ia - 1);
    G4String rems = aCommand.substr(ibx + 1, aCommand.length() - ibx);
    G4String* alVal = aliasList->FindAlias(alis);
    if(alVal == nullptr)
    {
      G4cerr << "Alias <" << alis << "> not found -- command ignored" << G4endl;
      G4String nullStr;
      return nullStr;
    }
    aCommand = subs + (*alVal) + rems;
    ia       = aCommand.find("{");
  }
  return aCommand;
}

// --------------------------------------------------------------------
G4int G4UImanager::ApplyCommand(const G4String& aCmd)
{
  return ApplyCommand(aCmd.data());
}

// --------------------------------------------------------------------
G4int G4UImanager::ApplyCommand(const char* aCmd)
{
  G4String aCommand = SolveAlias(aCmd);
  if(aCommand.empty())
  {
    return fAliasNotFound;
  }
  if(verboseLevel != 0)
  {
    if(isMaster)
    {
      fLastCommandOutputTreated = false;
    }
    G4cout << aCommand << G4endl;
  }
  G4String commandString;
  G4String commandParameter;

  std::size_t i = aCommand.find(" ");
  if(i != std::string::npos)
  {
    commandString    = aCommand.substr(0, i);
    commandParameter = aCommand.substr(i + 1, aCommand.length() - (i + 1));
  }
  else
  {
    commandString = aCommand;
  }

  // remove doubled slash
  std::size_t len = commandString.length();
  std::size_t ll  = 0;
  G4String a1;
  G4String a2;
  while(ll < len - 1)
  {
    if(commandString.substr(ll, 2) == "//")
    {
      if(ll == 0)
      {
        // Safe because index argument always 0
        commandString.erase(ll, 1);
      }
      else
      {
        a1            = commandString.substr(0, ll);
        a2            = commandString.substr(ll + 1, len - ll - 1);
        commandString = a1 + a2;
      }
      --len;
    }
    else
    {
      ++ll;
    }
  }

  if(isMaster && bridges != nullptr)
  {
    for(auto bridge : *bridges)
    {
      G4int leng = bridge->DirLength();
      if(commandString.substr(0, leng) == bridge->DirName())
      {
        return bridge->LocalUI()->ApplyCommand(commandString + " " +
                                               commandParameter);
      }
    }
  }

  G4UIcommand* targetCommand = treeTop->FindPath(commandString);
  if(targetCommand == nullptr)
  {
    if(ignoreCmdNotFound)
    {
      if(stackCommandsForBroadcast)
      {
        commandStack->push_back(commandString + " " + commandParameter);
      }
      return fCommandSucceeded;
    }
    else
    {
      return fCommandNotFound;
    }
  }

  if(stackCommandsForBroadcast && targetCommand->ToBeBroadcasted())
  {
    commandStack->push_back(commandString + " " + commandParameter);
  }

  if(!(targetCommand->IsAvailable()))
  {
    return fIllegalApplicationState;
  }

  if(saveHistory)
  {
    historyFile << aCommand << G4endl;
  }
  if(G4int(histVec.size()) >= maxHistSize)
  {
    histVec.erase(histVec.begin());
  }
  histVec.push_back(aCommand);

  targetCommand->ResetFailure();
  G4int commandFailureCode = targetCommand->DoIt(commandParameter);
  if(commandFailureCode == 0)
  {
    G4int additionalFailureCode = targetCommand->IfCommandFailed();
    if(additionalFailureCode > 0)
    {
      G4ExceptionDescription msg;
      msg << targetCommand->GetFailureDescription() << "\n"
          << "Error code : " << additionalFailureCode;
      G4Exception("G4UImanager::ApplyCommand", "UIMAN0123", JustWarning, msg);
      commandFailureCode += additionalFailureCode;
    }
  }
  return commandFailureCode;
}


// --------------------------------------------------------------------
G4UIcommand* G4UImanager::FindCommand(const G4String& aCmd)
{
  return FindCommand(aCmd.data());
}

// --------------------------------------------------------------------
G4UIcommand* G4UImanager::FindCommand(const char* aCmd)
{
  G4String aCommand = SolveAlias(aCmd);
  if(aCommand.empty()) 
  { return nullptr; }

  G4String commandString;

  std::size_t i = aCommand.find(" ");
  if(i != std::string::npos)
  { commandString = aCommand.substr(0, i); }
  else
  { commandString = aCommand; }

  return treeTop->FindPath(commandString);
}  

// --------------------------------------------------------------------
void G4UImanager::StoreHistory(const char* fileName)
{
  StoreHistory(true, fileName);
}

// --------------------------------------------------------------------
void G4UImanager::StoreHistory(G4bool historySwitch, const char* fileName)
{
  if(historySwitch)
  {
    if(saveHistory)
    {
      historyFile.close();
    }
    historyFile.open((char*) fileName);
    saveHistory = true;
  }
  else
  {
    historyFile.close();
    saveHistory = false;
  }
  saveHistory = historySwitch;
}

// --------------------------------------------------------------------
void G4UImanager::PauseSession(const char* msg)
{
  if(session != nullptr)
  {
    session->PauseSessionStart(msg);
  }
}

// --------------------------------------------------------------------
void G4UImanager::ListCommands(const char* direct)
{
  G4UIcommandTree* comTree = FindDirectory(direct);
  if(comTree != nullptr)
  {
    comTree->List();
  }
  else
  {
    G4cout << direct << " is not found." << G4endl;
  }
}

// --------------------------------------------------------------------
G4UIcommandTree* G4UImanager::FindDirectory(const char* dirName)
{
  G4String aDirName  = dirName;
  G4String targetDir = G4StrUtil::strip_copy(aDirName);
  if(targetDir.back() != '/')
  {
    targetDir += "/";
  }
  G4UIcommandTree* comTree = treeTop;
  if(targetDir == "/")
  {
    return comTree;
  }
  std::size_t idx = 1;
  while(idx < targetDir.length() - 1)
  {
    std::size_t i            = targetDir.find("/", idx);
    G4String targetDirString = targetDir.substr(0, i + 1);
    comTree                  = comTree->GetTree(targetDirString);
    if(comTree == nullptr)
    {
      return nullptr;
    }
    idx = i + 1;
  }
  return comTree;
}

// --------------------------------------------------------------------
G4bool G4UImanager::Notify(G4ApplicationState requestedState)
{
  if(pauseAtBeginOfEvent)
  {
    if(requestedState == G4State_EventProc &&
       G4StateManager::GetStateManager()->GetPreviousState() ==
         G4State_GeomClosed)
    {
      PauseSession("BeginOfEvent");
    }
  }
  if(pauseAtEndOfEvent)
  {
    if(requestedState == G4State_GeomClosed &&
       G4StateManager::GetStateManager()->GetPreviousState() ==
         G4State_EventProc)
    {
      PauseSession("EndOfEvent");
    }
  }
  return true;
}

// --------------------------------------------------------------------
void G4UImanager::SetCoutDestination(G4UIsession* const value)
{
  G4coutbuf.SetDestination(value);
  G4cerrbuf.SetDestination(value);
}

// --------------------------------------------------------------------
void G4UImanager::SetAlias(const char* aliasLine)
{
  G4String aLine      = aliasLine;
  std::size_t i       = aLine.find(" ");
  G4String aliasName  = aLine.substr(0, i);
  G4String aliasValue = aLine.substr(i + 1, aLine.length() - (i + 1));
  if(aliasValue[0] == '"')
  {
    G4String strippedValue;
    if(aliasValue.back() == '"')
    {
      strippedValue = aliasValue.substr(1, aliasValue.length() - 2);
    }
    else
    {
      strippedValue = aliasValue.substr(1, aliasValue.length() - 1);
    }
    aliasValue = strippedValue;
  }

  aliasList->ChangeAlias(aliasName, aliasValue);
}

// --------------------------------------------------------------------
void G4UImanager::RemoveAlias(const char* aliasName)
{
  G4String aL          = aliasName;
  G4String targetAlias = G4StrUtil::strip_copy(aL);
  aliasList->RemoveAlias(targetAlias);
}

// --------------------------------------------------------------------
void G4UImanager::ListAlias()
{
  aliasList->List();
}

// --------------------------------------------------------------------
void G4UImanager::CreateHTML(const char* dir)
{
  G4UIcommandTree* tr = FindDirectory(dir);
  if(tr != nullptr)
  {
    tr->CreateHTML();
  }
  else
  {
    G4cerr << "Directory <" << dir << "> is not found." << G4endl;
  }
}

// --------------------------------------------------------------------
void G4UImanager::ParseMacroSearchPath()
{
  searchDirs.clear();

  std::size_t idxfirst = 0;
  std::size_t idxend   = 0;
  G4String pathstring  = "";
  while((idxend = searchPath.find(':', idxfirst)) != G4String::npos)
  {
    pathstring = searchPath.substr(idxfirst, idxend - idxfirst);
    if(!pathstring.empty())
    {
      searchDirs.push_back(pathstring);
    }
    idxfirst = idxend + 1;
  }

  pathstring = searchPath.substr(idxfirst, searchPath.size() - idxfirst);
  if(!pathstring.empty())
  {
    searchDirs.push_back(pathstring);
  }
}

// --------------------------------------------------------------------
static G4bool FileFound(const G4String& fname)
{
  G4bool qopen = false;
  std::ifstream fs;
  fs.open(fname.c_str(), std::ios::in);
  if(fs.good())
  {
    fs.close();
    qopen = true;
  }
  return qopen;
}

// --------------------------------------------------------------------
G4String G4UImanager::FindMacroPath(const G4String& fname) const
{
  G4String macrofile = fname;

  for(const auto& searchDir : searchDirs)
  {
    G4String fullpath = searchDir + "/" + fname;
    if(FileFound(fullpath))
    {
      macrofile = fullpath;
      break;
    }
  }
  return macrofile;
}

// --------------------------------------------------------------------
std::vector<G4String>* G4UImanager::GetCommandStack()
{
  std::vector<G4String>* returnValue = commandStack;
  commandStack                       = new std::vector<G4String>;
  return returnValue;
}

// --------------------------------------------------------------------
void G4UImanager::RegisterBridge(G4UIbridge* brg)
{
  if(brg->LocalUI() == this)
  {
    G4Exception("G4UImanager::RegisterBridge()", "UI7002", FatalException,
                "G4UIBridge cannot bridge between same object.");
  }
  else
  {
    bridges->push_back(brg);
  }
}

// --------------------------------------------------------------------
void G4UImanager::SetUpForAThread(G4int tId)
{
  threadID = tId;
  G4iosInitialization();
  threadCout = new G4MTcoutDestination(threadID);
  threadCout->SetIgnoreCout(igThreadID);
}

// --------------------------------------------------------------------
void G4UImanager::SetUpForSpecialThread(const G4String& pref)
{
  threadID = G4Threading::GENERICTHREAD_ID;
  G4Threading::G4SetThreadId(threadID);
  G4iosInitialization();
  threadCout = new G4MTcoutDestination(threadID);
  threadCout->SetPrefixString(pref);
  threadCout->SetIgnoreCout(igThreadID);
}

// --------------------------------------------------------------------
void G4UImanager::SetCoutFileName(const G4String& fileN, G4bool ifAppend)
{
  // for sequential mode, ignore this method.
  if(threadID < 0)
  {
    return;
  }

  if(fileN == "**Screen**")
  {
    threadCout->SetCoutFileName(fileN, ifAppend);
  }
  else
  {
    std::stringstream fn;
    fn << "G4W_" << threadID << "_" << fileN;
    threadCout->SetCoutFileName(fn.str(), ifAppend);
  }
}

// --------------------------------------------------------------------
void G4UImanager::SetCerrFileName(const G4String& fileN, G4bool ifAppend)
{
  // for sequential mode, ignore this method.
  if(threadID < 0)
  {
    return;
  }

  if(fileN == "**Screen**")
  {
    threadCout->SetCerrFileName(fileN, ifAppend);
  }
  else
  {
    std::stringstream fn;
    fn << "G4W_" << threadID << "_" << fileN;
    threadCout->SetCerrFileName(fn.str(), ifAppend);
  }
}

// --------------------------------------------------------------------
void G4UImanager::SetThreadPrefixString(const G4String& s)
{
  // for sequential mode, ignore this method.
  if(threadID < 0)
  {
    return;
  }
  threadCout->SetPrefixString(s);
}

// --------------------------------------------------------------------
void G4UImanager::SetThreadUseBuffer(G4bool flg)
{
  // for sequential mode, ignore this method.
  if(threadID < 0)
  {
    return;
  }
  threadCout->EnableBuffering(flg);
}

// --------------------------------------------------------------------
void G4UImanager::SetThreadIgnore(G4int tid)
{
  // for sequential mode, ignore this method.
  if(threadID < 0)
  {
    igThreadID = tid;
    return;
  }
  threadCout->SetIgnoreCout(tid);
}

// --------------------------------------------------------------------
void G4UImanager::SetThreadIgnoreInit(G4bool flg)
{
  // for sequential mode, ignore this method.
  if(threadID < 0)
  {
    return;
  }
  threadCout->SetIgnoreInit(flg);
}

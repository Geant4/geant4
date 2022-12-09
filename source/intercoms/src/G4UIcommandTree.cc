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
// G4UIcommandTree
//
// Author: Makoto Asai (SLAC), 1998
// Midified: Makoto Asai (SLAC), 2021
//   Improve output HTML file layout and add option to sort
//   command/directory names in alphabetic order
// --------------------------------------------------------------------

#include "G4UIcommandTree.hh"
#include "G4UIdirectory.hh"
#include "G4StateManager.hh"
#include "G4UImanager.hh"
#include <fstream>
#include "G4ios.hh"

// --------------------------------------------------------------------
G4UIcommandTree::G4UIcommandTree(const char* thePathName)
{
  pathName = thePathName;
}

// --------------------------------------------------------------------
G4UIcommandTree::~G4UIcommandTree()
{
  for(std::size_t i = 0; i < tree.size(); ++i)
  {
    delete tree[i];
  }
}

// --------------------------------------------------------------------
G4bool G4UIcommandTree::operator==(const G4UIcommandTree& right) const
{
  return (pathName == right.GetPathName());
}

// --------------------------------------------------------------------
G4bool G4UIcommandTree::operator!=(const G4UIcommandTree& right) const
{
  return (pathName != right.GetPathName());
}

// --------------------------------------------------------------------
void G4UIcommandTree::AddNewCommand(G4UIcommand* newCommand,
                                    G4bool workerThreadOnly)
{
  G4String commandPath   = newCommand->GetCommandPath();
  G4String remainingPath = commandPath;
  remainingPath.erase(0, pathName.length());
  if(remainingPath.empty())
  {
    if(guidance == nullptr)
    {
      guidance = newCommand;
      if(!(newCommand->ToBeBroadcasted()))
      { broadcastCommands = false; }
      if(workerThreadOnly)
      { newCommand->SetWorkerThreadOnly(); }
    }
    return;
  }

  if(guidance != nullptr)
  {
    auto* dir = static_cast<G4UIdirectory*>(guidance);
    ifSort = dir->IfSort();
  }
  std::size_t i = remainingPath.find('/');
  if(i == std::string::npos)
  {
    // Adding a new command to this directory
    std::size_t n_commandEntry = command.size();
    for(std::size_t i_thCommand = 0; i_thCommand < n_commandEntry; ++i_thCommand)
    {
      if(remainingPath == command[i_thCommand]->GetCommandName())
      {
        // a command of same name has already defined. do nothing and return.
        if(G4UImanager::GetUIpointer()->GetVerboseLevel() > 8)
        {
          G4ExceptionDescription ed;
          ed << "Command <" << commandPath << "> already exist. New command is not added.";
          G4Exception("G4UIcommandTree::AddNewCommand","UI_ComTree_001",
                     //FatalException,
                     JustWarning,
                     ed);
        }
        return;
      }
    }
    if(!broadcastCommands)
    { newCommand->SetToBeBroadcasted(false); }
    if(workerThreadOnly)
    { newCommand->SetWorkerThreadOnly(); }
    if(ifSort)
    {
      auto j = command.cbegin();
      for(; j != command.cend(); ++j) {
        if (newCommand->GetCommandPath() < (*j)->GetCommandPath()) { break; }
      }
      command.insert(j,newCommand);
    }
    else
    { command.push_back(newCommand); }
    return;
  }
  else
  {
    // Adding a new command to a sub-directory
    G4String nextPath = pathName;
    nextPath.append(remainingPath.substr(0, i + 1));
    std::size_t n_treeEntry = tree.size();
    for(std::size_t i_thTree = 0; i_thTree < n_treeEntry; ++i_thTree)
    {
      if(nextPath == tree[i_thTree]->GetPathName())
      {
        if(!broadcastCommands)
        {
          newCommand->SetToBeBroadcasted(false);
        }
        tree[i_thTree]->AddNewCommand(newCommand, workerThreadOnly);
        return;
      }
    }
    // Creating a new sub-directory
    auto* newTree = new G4UIcommandTree(nextPath);
    if(ifSort)
    {
      auto j = tree.cbegin();
      for(; j != tree.cend(); ++j) {
        if (newTree->GetPathName() < (*j)->GetPathName()) { break; }
      }
      tree.insert(j,newTree);
    }
    else
    { tree.push_back(newTree); }
    if(!broadcastCommands)
    { newCommand->SetToBeBroadcasted(false); }
    // In case a new sub-directry is created with a new G4UIdirectory
    // (most-likely this is the case), inherit the sort flag
    newCommand->SetDefaultSortFlag(ifSort);
    newTree->AddNewCommand(newCommand, workerThreadOnly);
    return;
  }
}

// --------------------------------------------------------------------
void G4UIcommandTree::RemoveCommand(G4UIcommand* aCommand,
                                    G4bool workerThreadOnly)
{
  if(workerThreadOnly && !(aCommand->IsWorkerThreadOnly()))
  {
    return;
  }
  G4String commandPath   = aCommand->GetCommandPath();
  G4String remainingPath = commandPath;
  remainingPath.erase(0, pathName.length());
  if(remainingPath.empty())
  {
    guidance = nullptr;
  }
  else
  {
    std::size_t i = remainingPath.find('/');
    if(i == std::string::npos)
    {
      // Find command
      std::size_t n_commandEntry = command.size();
      for(std::size_t i_thCommand = 0; i_thCommand < n_commandEntry; ++i_thCommand)
      {
        if(remainingPath == command[i_thCommand]->GetCommandName())
        {
          command.erase(command.begin() + i_thCommand);
          break;
        }
      }
    }
    else
    {
      // Find path
      G4String nextPath = pathName;
      nextPath.append(remainingPath.substr(0, i + 1));
      std::size_t n_treeEntry = tree.size();
      for(std::size_t i_thTree = 0; i_thTree < n_treeEntry; ++i_thTree)
      {
        if(nextPath == tree[i_thTree]->GetPathName())
        {
          tree[i_thTree]->RemoveCommand(aCommand);
          G4int n_commandRemain = tree[i_thTree]->GetCommandEntry();
          G4int n_treeRemain    = tree[i_thTree]->GetTreeEntry();
          if(n_commandRemain == 0 && n_treeRemain == 0)
          {
            G4UIcommandTree* emptyTree = tree[i_thTree];
            tree.erase(tree.begin() + i_thTree);
            delete emptyTree;
          }
          break;
        }
      }
    }
  }
}

// --------------------------------------------------------------------
G4UIcommand* G4UIcommandTree::FindPath(const char* commandPath) const
{
  // This function tries to match a command name

  G4String remainingPath = commandPath;
  if(remainingPath.find(pathName) == std::string::npos)
  {
    return nullptr;
  }
  remainingPath.erase(0, pathName.length());
  std::size_t i = remainingPath.find('/');
  if(i == std::string::npos)
  {
    // Find command
    std::size_t n_commandEntry = command.size();
    for(std::size_t i_thCommand = 0; i_thCommand < n_commandEntry; ++i_thCommand)
    {
      if(remainingPath == command[i_thCommand]->GetCommandName())
      {
        return command[i_thCommand];
      }
    }
  }
  else
  {
    // Find path
    G4String nextPath = pathName;
    nextPath.append(remainingPath.substr(0, i + 1));
    std::size_t n_treeEntry = tree.size();
    for(std::size_t i_thTree = 0; i_thTree < n_treeEntry; ++i_thTree)
    {
      if(nextPath == tree[i_thTree]->GetPathName())
      {
        return tree[i_thTree]->FindPath(commandPath);
      }
    }
  }
  return nullptr;
}

// --------------------------------------------------------------------
G4UIcommandTree* G4UIcommandTree::FindCommandTree(const char* commandPath)
{
  // Try to match a command or a path with the one given.
  // @commandPath : command or path to match
  // @return the commandTree found or nullptr if not

  G4String remainingPath = commandPath;
  if(remainingPath.find(pathName) == std::string::npos)
  {
    return nullptr;
  }
  remainingPath.erase(0, pathName.length());
  std::size_t i = remainingPath.find('/');
  if(i != std::string::npos)
  {
    // Find path
    G4String nextPath = pathName;
    nextPath.append(remainingPath.substr(0, i + 1));
    std::size_t n_treeEntry = tree.size();
    for(std::size_t i_thTree = 0; i_thTree < n_treeEntry; ++i_thTree)
    {
      if(tree[i_thTree]->GetPathName() == commandPath)
      {
        return tree[i_thTree];
      }
      else if(nextPath == tree[i_thTree]->GetPathName())
      {
        return tree[i_thTree]->FindCommandTree(commandPath);
      }
    }
  }
  else
  {
    return this;
  }
  return nullptr;
}

// --------------------------------------------------------------------
G4String G4UIcommandTree::CompleteCommandPath(const G4String& aCommandPath)
{
  G4String pName         = aCommandPath;
  G4String remainingPath = aCommandPath;
  G4String empty         = "";
  G4String matchingPath  = empty;

  // find the tree
  auto jpre = pName.rfind('/');
  if(jpre != G4String::npos)
  {
    pName.erase(jpre + 1);
  }
  G4UIcommandTree* aTree = FindCommandTree(pName);

  if(aTree == nullptr)
  {
    return empty;
  }

  if(pName.find(pName) == std::string::npos)
  {
    return empty;
  }

  std::vector<G4String> paths;

  // list matched directories/commands
  G4String strtmp;
  G4int nMatch = 0;

  G4int Ndir = aTree->GetTreeEntry();
  G4int Ncmd = aTree->GetCommandEntry();

  // directory ...
  for(G4int idir = 1; idir <= Ndir; ++idir)
  {
    G4String fpdir = aTree->GetTree(idir)->GetPathName();
    // matching test
    if(fpdir.find(remainingPath, 0) == 0)
    {
      if(nMatch == 0)
      {
        matchingPath = fpdir;
      }
      else
      {
        matchingPath = GetFirstMatchedString(fpdir, matchingPath);
      }
      ++nMatch;
      paths.push_back(fpdir);
    }
  }

  if(paths.size() >= 2)
  {
    G4cout << "Matching directories :" << G4endl;
    for(const auto& path : paths)
    {
      G4cout << path << G4endl;
    }
  }

  // command ...
  std::vector<G4String> commands;

  for(G4int icmd = 1; icmd <= Ncmd; ++icmd)
  {
    G4String fpcmd =
      aTree->GetPathName() + aTree->GetCommand(icmd)->GetCommandName();
    // matching test
    if(fpcmd.find(remainingPath, 0) == 0)
    {
      if(nMatch == 0)
      {
        matchingPath = fpcmd + " ";
      }
      else
      {
        strtmp       = fpcmd + " ";
        matchingPath = GetFirstMatchedString(matchingPath, strtmp);
      }
      nMatch++;
      commands.emplace_back(fpcmd + " ");
    }
  }

  if(commands.size() >= 2)
  {
    G4cout << "Matching commands :" << G4endl;
    for(const auto& matched : commands)
    {
      G4cout << matched << G4endl;
    }
  }

  return matchingPath;
}

// --------------------------------------------------------------------
G4String G4UIcommandTree::GetFirstMatchedString(const G4String& str1,
                                                const G4String& str2) const
{
  std::size_t nlen1 = str1.length();
  std::size_t nlen2 = str2.length();

  std::size_t nmin = nlen1 < nlen2 ? nlen1 : nlen2;

  G4String strMatched;
  for(G4int i = 0; i < (G4int)nmin; ++i)
  {
    if(str1[i] == str2[i])
    {
      strMatched += str1[i];
    }
    else
    {
      break;
    }
  }

  return strMatched;
}

// --------------------------------------------------------------------
void G4UIcommandTree::ListCurrent() const
{
  G4cout << "Command directory path : " << pathName << G4endl;
  if(guidance != nullptr)
  {
    guidance->List();
  }
  G4cout << " Sub-directories : " << G4endl;
  std::size_t n_treeEntry = tree.size();
  for(std::size_t i_thTree = 0; i_thTree < n_treeEntry; ++i_thTree)
  {
    G4cout << "   " << tree[i_thTree]->GetPathName();
    if((tree[i_thTree]->GetGuidance() != nullptr) &&
       tree[i_thTree]->GetGuidance()->IsWorkerThreadOnly())
    {
      G4cout << " @ ";
    }
    else
    {
      G4cout << "   ";
    }
    G4cout << tree[i_thTree]->GetTitle() << G4endl;
  }
  G4cout << " Commands : " << G4endl;
  std::size_t n_commandEntry = command.size();
  for(std::size_t i_thCommand = 0; i_thCommand < n_commandEntry; ++i_thCommand)
  {
    G4cout << "   " << command[i_thCommand]->GetCommandName();
    if(command[i_thCommand]->IsWorkerThreadOnly())
    {
      G4cout << " @ ";
    }
    else
    {
      G4cout << " * ";
    }
    G4cout << command[i_thCommand]->GetTitle() << G4endl;
  }
}

// --------------------------------------------------------------------
void G4UIcommandTree::ListCurrentWithNum() const
{
  G4cout << "Command directory path : " << pathName << G4endl;
  if(guidance != nullptr)
  {
    guidance->List();
  }
  G4int i = 0;
  G4cout << " Sub-directories : " << G4endl;
  std::size_t n_treeEntry = tree.size();
  for(std::size_t i_thTree = 0; i_thTree < n_treeEntry; ++i_thTree)
  {
    ++i;
    G4cout << " " << i << ") " << tree[i_thTree]->GetPathName() << "   "
           << tree[i_thTree]->GetTitle() << G4endl;
  }
  G4cout << " Commands : " << G4endl;
  std::size_t n_commandEntry = command.size();
  for(std::size_t i_thCommand = 0; i_thCommand < n_commandEntry; ++i_thCommand)
  {
    ++i;
    G4cout << " " << i << ") " << command[i_thCommand]->GetCommandName()
           << " * " << command[i_thCommand]->GetTitle() << G4endl;
  }
}

// --------------------------------------------------------------------
void G4UIcommandTree::List() const
{
  ListCurrent();
  std::size_t n_commandEntry = command.size();
  for(std::size_t i_thCommand = 0; i_thCommand < n_commandEntry; ++i_thCommand)
  {
    command[i_thCommand]->List();
  }
  std::size_t n_treeEntry = tree.size();
  for(std::size_t i_thTree = 0; i_thTree < n_treeEntry; ++i_thTree)
  {
    tree[i_thTree]->List();
  }
}

// --------------------------------------------------------------------
G4String G4UIcommandTree::CreateFileName(const char* pName)
{
  G4String fn = pName;
  std::size_t idxs;
  while((idxs = fn.find("/")) != std::string::npos)
  {
    fn[(G4int)idxs] = '_';
  }
  fn += ".html";
  return fn;
}

// --------------------------------------------------------------------
G4String G4UIcommandTree::ModStr(const char* strS)
{
  G4String sx;
  G4String str = strS;
  for(G4int i = 0; i < G4int(str.length()); ++i)
  {
    char c = str[i];
    switch(c)
    {
      case '<':
        sx += "&lt;";
        break;
      case '>':
        sx += "&gt;";
        break;
      case '&':
        sx += "&amp;";
        break;
      default:
        sx += c;
    }
  }
  return sx;
}

// --------------------------------------------------------------------
void G4UIcommandTree::CreateHTML(const G4String& sideBar)
{
  G4String ofileName = CreateFileName(pathName);
  std::ofstream oF(ofileName, std::ios::out);

  oF << "<html><head><title>Commands in " << ModStr(pathName)
     << "</title></head>" << G4endl;
    oF << "<style> \
    table,table td,table th { \
       border:1px solid #eee \
    } \
    table td,table th { \
      padding:5px 20px; \
      line-height:1.3; \
      text-align:inherit \
    } \
    a { \
      color:#17a81a; \
      text-decoration:none; \
      transition-duration:0.3s \
    } \
      a:hover { \
      color:#17a81a \
    } \
    table { \
      border-collapse:collapse; \
      border-spacing:0; \
      margin-bottom:5px; \
    } \
    h1 { \
      font-size:2.25em; \
      font-weight:300; \
      letter-spacing:-1px; \
      line-height:1.15em; \
      margin-bottom:0.5em; \
      word-wrap:break-word \
    } \
    h2 { \
      font-size:1.5em; \
      font-weight:300; \
      letter-spacing:-1px; \
      line-height:1.15em; \
      margin-bottom:0.5em; \
      word-wrap:break-word \
    } \
    h3 { \
      color:#26282a; \
      font-weight:300; \
      font-size:1.3em; \
      padding:15px 0 15px 0; \
      border-bottom:2px #eee solid; \
      word-wrap:break-word \
    } \
    .sidebar { \
      display:block; \
      position:relative; \
      position:sticky; \
      float:left; \
      -webkit-box-sizing:border-box; \
      -moz-box-sizing:border-box; \
      -ms-box-sizing:border-box; \
      box-sizing:border-box; \
      width:20%; \
      padding-right:20px \
    } \
    .context { \
    width:80%; \
    display:inline-block; \
    background-color:#fff; \
    padding: 25px 35px 20px 30px; \
    -webkit-box-sizing:border-box; \
    -moz-box-sizing:border-box; \
    -ms-box-sizing:border-box; \
    box-sizing:border-box \
  } \
    </style>"<< G4endl;
  oF << "<body bgcolor=\"#ffffff\">" << G4endl;

    // Left Panel
  if (createHTMLTreeLevel == 0 ) {
    oF << "<div class=\"sidebar\">" << sideBar << "</div>" << G4endl;
  }
  // Right Panel
  oF << "<div class=\"context\">";
  oF <<  "<h1>" << ModStr(pathName) << "</h1>" << G4endl;

  if(guidance != nullptr)
  {
    for(G4int i = 0; i < (G4int)guidance->GetGuidanceEntries(); ++i)
    {
      oF << ModStr(guidance->GetGuidanceLine(i)) << "<br>" << G4endl;
    }
  }
  if(!tree.empty())
  {
    G4String menu = "";
    G4String newSideBar = "";
    menu += "<h2>Sub-directories </h2><table>";
      newSideBar += "<h2><a href=\"" + ofileName + "\">Top level </a></h2><table>";
    // Build menu short version
      for(auto& i_thTree : tree)
      {
        newSideBar += "<tr><td><a href=\"" +
                      CreateFileName(i_thTree->GetPathName()) + "\">" +
                      ModStr(i_thTree->GetPathName()) + "</a>";
    }
    // Build menu
    for(auto& i_thTree : tree)
    {
      menu += "<tr><td><a href=\"" + CreateFileName(i_thTree->GetPathName()) +
              "\">" + ModStr(i_thTree->GetPathName()) + "</a>";
      menu += "</td><td>" + ModStr(i_thTree->GetTitle()) + "</tr>";
    }
    menu += "</table>";
      newSideBar += "</table>";
      for(auto& i_thTree : tree)
      {
        createHTMLTreeLevel ++;
        i_thTree->CreateHTML(newSideBar);
        createHTMLTreeLevel --;
    }
    oF << menu << G4endl;
  }

  if(!command.empty())
  {
    oF << "<h2>Commands </h2>" << G4endl;
      
    // resume
    oF << "<table>" << G4endl;
    for(std::size_t i_thCommand = 0; i_thCommand < command.size(); ++i_thCommand)
    {
      G4UIcommand* cmd = command[i_thCommand];
      oF << "<tr><td><a href=\"#c"<< i_thCommand << "\">"<< ModStr(cmd->GetCommandName());
      oF << "</a></td></tr>" << G4endl;
    }
    oF << "</table>" << G4endl;
    for(std::size_t i_thCommand = 0; i_thCommand < command.size(); ++i_thCommand)
    {
      G4UIcommand* cmd = command[i_thCommand];
      oF << "<h3 id=\"c" << i_thCommand << "\">" << ModStr(cmd->GetCommandName());
      if(cmd->GetParameterEntries() > 0)
      {
        for(G4int i_thParam = 0; i_thParam < (G4int)cmd->GetParameterEntries();
            ++i_thParam)
        {
          oF << " [<i>"
             << ModStr(cmd->GetParameter(i_thParam)->GetParameterName())
             << "</i>]";
        }
      }
      oF << "</h3>" << G4endl;
      oF << "<p>" << G4endl;
      for(G4int i = 0; i < (G4int)cmd->GetGuidanceEntries(); ++i)
      {
        oF << ModStr(cmd->GetGuidanceLine(i)) << "<br>" << G4endl;
      }
      if(!(cmd->GetRange()).empty())
      {
        oF << "<p>Range : " << ModStr(cmd->GetRange()) << G4endl;
      }
      std::vector<G4ApplicationState>* availabelStateList = cmd->GetStateList();
      if(availabelStateList->size() == 6)
      {
        oF << "<p>Available at all Geant4 states." << G4endl;
      }
      else
      {
        oF << "<p>Available Geant4 state(s) : ";
        for(auto& ias : *availabelStateList)
        {
          oF << G4StateManager::GetStateManager()->GetStateString(ias) << " "
             << G4endl;
        }
      }
      if(cmd->GetParameterEntries() > 0)
      {
        oF << "<p>Parameters<table border=1>" << G4endl;
        for(G4int i_thParam = 0; i_thParam < (G4int)cmd->GetParameterEntries();
            ++i_thParam)
        {
          G4UIparameter* prm = cmd->GetParameter(i_thParam);
          oF << "<tr><td>" << ModStr(prm->GetParameterName()) << G4endl;
          oF << "<td>type " << prm->GetParameterType() << G4endl;
          oF << "<td>";
          if(prm->IsOmittable())
          {
            oF << "Omittable : ";
            if(prm->GetCurrentAsDefault())
            {
              oF << "current value is used as the default value." << G4endl;
            }
            else
            {
              oF << "default value = " << prm->GetDefaultValue() << G4endl;
            }
          }
          oF << "<td>";
          if(!(prm->GetParameterRange()).empty())
          {
            oF << "Parameter range : " << ModStr(prm->GetParameterRange())
               << G4endl;
          }
          else if(!(prm->GetParameterCandidates()).empty())
          {
            oF << "Parameter candidates : "
               << ModStr(prm->GetParameterCandidates()) << G4endl;
          }
        }
        oF << "</table>" << G4endl;
      }
    }
  }
  oF << "</div></body></html>" << G4endl;
  oF.close();
}

// --------------------------------------------------------------------
G4UIcommandTree* G4UIcommandTree::GetTree(const char* comNameC)
{
  G4String comName = comNameC;
  for(auto& i : tree)
  {
    if(comName == i->GetPathName())
    {
      return i;
    }
  }
  return nullptr;
}

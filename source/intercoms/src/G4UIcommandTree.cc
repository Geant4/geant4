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
//
// $Id: G4UIcommandTree.cc 77651 2013-11-27 08:47:55Z gcosmo $
//

#include "G4UIcommandTree.hh"
#include "G4StateManager.hh"
#include <fstream>
#include "G4ios.hh"

G4UIcommandTree::G4UIcommandTree()
:guidance(NULL),broadcastCommands(true)
{ }

G4UIcommandTree::G4UIcommandTree(const char * thePathName)
:guidance(NULL),broadcastCommands(true)
{
  pathName = thePathName;
}

G4UIcommandTree::~G4UIcommandTree()
{
  G4int i;
  G4int n_treeEntry = tree.size();
  for( i=0; i < n_treeEntry; i++ )
  { delete tree[i]; }
}

G4int G4UIcommandTree::operator==(const G4UIcommandTree &right) const
{
  return ( pathName == right.GetPathName() );
}

G4int G4UIcommandTree::operator!=(const G4UIcommandTree &right) const
{
  return ( pathName != right.GetPathName() );
}

void G4UIcommandTree::AddNewCommand(G4UIcommand *newCommand, G4bool workerThreadOnly)
{
  G4String commandPath = newCommand->GetCommandPath();
  G4String remainingPath = commandPath;
  remainingPath.remove(0,pathName.length());
  if( remainingPath.isNull() )
  {
    if(!guidance)
    {
      guidance = newCommand;
      if(!(newCommand->ToBeBroadcasted())) broadcastCommands = false;
      if(workerThreadOnly) newCommand->SetWorkerThreadOnly();
    }
    return;
  }
  G4int i = remainingPath.first('/');
  if( i == G4int(std::string::npos) )
  {
    // Find command
    G4int n_commandEntry = command.size();
    for( G4int i_thCommand = 0; i_thCommand < n_commandEntry; i_thCommand++ )
    {
      if( remainingPath == command[i_thCommand]->GetCommandName() )
      { return; }
    }
    if(!broadcastCommands) newCommand->SetToBeBroadcasted(false);
    if(workerThreadOnly) newCommand->SetWorkerThreadOnly();
    command.push_back( newCommand );
    return;
  }
  else
  {
    // Find path
    G4String nextPath = pathName;
    nextPath.append(remainingPath(0,i+1));
    G4int n_treeEntry = tree.size();
    for( G4int i_thTree = 0; i_thTree < n_treeEntry; i_thTree++ )
    {
      if( nextPath == tree[i_thTree]->GetPathName() )
      { 
        if(!broadcastCommands) newCommand->SetToBeBroadcasted(false);
	tree[i_thTree]->AddNewCommand( newCommand, workerThreadOnly );
	return; 
      }
    }
    G4UIcommandTree * newTree = new G4UIcommandTree( nextPath );
    tree.push_back( newTree );
    if(!broadcastCommands) newCommand->SetToBeBroadcasted(false);
    newTree->AddNewCommand( newCommand, workerThreadOnly );
    return;
  }
}

void G4UIcommandTree::RemoveCommand(G4UIcommand *aCommand, G4bool workerThreadOnly)
{
  if(workerThreadOnly && !(aCommand->IsWorkerThreadOnly())) return;
  G4String commandPath = aCommand->GetCommandPath();
  G4String remainingPath = commandPath;
  remainingPath.remove(0,pathName.length());
  if( remainingPath.isNull() )
  {
    guidance = NULL;
  }
  else
  {
    G4int i = remainingPath.first('/');
    if( i == G4int(std::string::npos) )
    {
      // Find command
      G4int n_commandEntry = command.size();
      for( G4int i_thCommand = 0; i_thCommand < n_commandEntry; i_thCommand++ )
      {
        if( remainingPath == command[i_thCommand]->GetCommandName() )
        { 
          command.erase(command.begin()+i_thCommand);
          break;
        }
      }
    }
    else
    {
      // Find path
      G4String nextPath = pathName;
      nextPath.append(remainingPath(0,i+1));
      G4int n_treeEntry = tree.size();
      for( G4int i_thTree = 0; i_thTree < n_treeEntry; i_thTree++ )
      {
        if( nextPath == tree[i_thTree]->GetPathName() )
        { 
    	  tree[i_thTree]->RemoveCommand( aCommand );
    	  G4int n_commandRemain = tree[i_thTree]->GetCommandEntry();
        G4int n_treeRemain = tree[i_thTree]-> GetTreeEntry();
    	  if(n_commandRemain == 0 && n_treeRemain == 0)
    	  {
    	    G4UIcommandTree * emptyTree = tree[i_thTree];
    	    tree.erase(tree.begin()+i_thTree);
    	    delete emptyTree;
    	  }
    	  break;
        }
      }
    }
  }
}

// L. Garnier 01.28.08 This function has not a good name. In fact, it try
// to match a command name, not a path. It should be rename as FindCommandName

G4UIcommand * G4UIcommandTree::FindPath(const char* commandPath) const
{
  G4String remainingPath = commandPath;
  if( remainingPath.index( pathName ) == std::string::npos )
  { return NULL; }
  remainingPath.remove(0,pathName.length());
  G4int i = remainingPath.first('/');
  if( i == G4int(std::string::npos) )
  {
    // Find command
    G4int n_commandEntry = command.size();
    for( G4int i_thCommand = 0; i_thCommand < n_commandEntry; i_thCommand++ )
    {
      if( remainingPath == command[i_thCommand]->GetCommandName() )
      { return command[i_thCommand]; }
    }
  }
  else
  {
    // Find path
    G4String nextPath = pathName;
    nextPath.append(remainingPath(0,i+1));
    G4int n_treeEntry = tree.size();
    for( G4int i_thTree = 0; i_thTree < n_treeEntry; i_thTree++ )
    {
      if( nextPath == tree[i_thTree]->GetPathName() )
      { return tree[i_thTree]->FindPath( commandPath ); }
    }
  }
  return NULL;
}


/**
 * Try to match a command or a path with the one given.
 * @commandPath : command or path to match
 * @return the commandTree found or NULL if not
 */
G4UIcommandTree* G4UIcommandTree::FindCommandTree(const char* commandPath)
{
  G4String remainingPath = commandPath;
  if( remainingPath.index( pathName ) == std::string::npos )
  { return NULL; }
  remainingPath.remove(0,pathName.length());
  G4int i = remainingPath.first('/');
  if( i != G4int(std::string::npos) )
  {
    // Find path
    G4String nextPath = pathName;
    nextPath.append(remainingPath(0,i+1));
    G4int n_treeEntry = tree.size();
    for( G4int i_thTree = 0; i_thTree < n_treeEntry; i_thTree++ )
    {
      if (tree[i_thTree]->GetPathName() == commandPath) {
        return tree[i_thTree];
      }
      else if( nextPath == tree[i_thTree]->GetPathName() ) {
        return tree[i_thTree]->FindCommandTree( commandPath );
      }
    }
  } else {
    return this;
  }
  return NULL;
}

G4String G4UIcommandTree::CompleteCommandPath(const G4String& aCommandPath)
{
  G4String pName = aCommandPath;
  G4String remainingPath = aCommandPath;
  G4String empty = "";
  G4String matchingPath = empty;

  // find the tree
  G4int jpre= pName.last('/');
  if(jpre != G4int(G4String::npos)) pName.remove(jpre+1);
  G4UIcommandTree* aTree = FindCommandTree(pName);

  if (!aTree) {
    return empty;
  }
  
  if( pName.index( pName ) == std::string::npos ) return empty;

  std::vector<G4String> paths;

  // list matched directories/commands
  G4String strtmp;
  G4int nMatch= 0;

  int Ndir= aTree-> GetTreeEntry();
  int Ncmd= aTree-> GetCommandEntry();
  
  // directory ...
  for(G4int idir=1; idir<=Ndir; idir++) {
    G4String fpdir= aTree-> GetTree(idir)-> GetPathName();
    // matching test
    if( fpdir.index(remainingPath, 0) == 0) {
      if(nMatch==0) {
        matchingPath = fpdir;
      } else {
        matchingPath = GetFirstMatchedString(fpdir,matchingPath);
      }
      nMatch++;
      paths.push_back(fpdir);
    }
  }
  
  if (paths.size()>=2) {
    G4cout << "Matching directories :" << G4endl; 
    for( unsigned int i_thCommand = 0; i_thCommand < paths.size(); i_thCommand++ ) {
      G4cout << paths[i_thCommand] << G4endl; 
    }
  }
  
  // command ...
  std::vector<G4String> commands;

  for(G4int icmd=1; icmd<=Ncmd; icmd++){
    G4String fpcmd= aTree-> GetPathName() +
                    aTree-> GetCommand(icmd) -> GetCommandName();
    // matching test
    if( fpcmd.index(remainingPath, 0) ==0) {
      if(nMatch==0) {
        matchingPath= fpcmd + " ";
      } else {
        strtmp= fpcmd + " ";
        matchingPath= GetFirstMatchedString(matchingPath, strtmp);
      }
      nMatch++;
      commands.push_back(fpcmd+" ");
    }
  }

  if (commands.size()>=2) {
    G4cout << "Matching commands :" << G4endl; 
    for( unsigned int i_thCommand = 0; i_thCommand < commands.size(); i_thCommand++ ) {
      G4cout << commands[i_thCommand] << G4endl; 
    }
  }

  return matchingPath;
}


////////////////////////////////////////////////////////////////////
G4String G4UIcommandTree::GetFirstMatchedString(const G4String& str1, 
					 const G4String& str2) const
////////////////////////////////////////////////////////////////////
{
  int nlen1= str1.length();
  int nlen2= str2.length();

  int nmin = nlen1<nlen2 ? nlen1 : nlen2;

  G4String strMatched;
  for(size_t i=0; G4int(i)<nmin; i++){
    if(str1[i]==str2[i]) {
      strMatched+= str1[i];
    } else {
      break;
    }
  }

  return strMatched;
}

void G4UIcommandTree::ListCurrent() const
{
  G4cout << "Command directory path : " << pathName << G4endl;
  if( guidance != NULL ) guidance->List();
  G4cout << " Sub-directories : " << G4endl;
  G4int n_treeEntry = tree.size();
  for( G4int i_thTree = 0; i_thTree < n_treeEntry; i_thTree++ )
  { 
    G4cout << "   " << tree[i_thTree]->GetPathName();
    if(tree[i_thTree]->GetGuidance() &&
       tree[i_thTree]->GetGuidance()->IsWorkerThreadOnly())
    { G4cout << " @ "; }
    else
    { G4cout << "   "; }
    G4cout << tree[i_thTree]->GetTitle() << G4endl;
  }
  G4cout << " Commands : " << G4endl;
  G4int n_commandEntry = command.size();
  for( G4int i_thCommand = 0; i_thCommand < n_commandEntry; i_thCommand++ )
  {
    G4cout << "   " << command[i_thCommand]->GetCommandName();
    if(command[i_thCommand]->IsWorkerThreadOnly())
    { G4cout << " @ "; }
    else
    { G4cout << " * "; }
    G4cout << command[i_thCommand]->GetTitle() << G4endl;
  }
}

void G4UIcommandTree::ListCurrentWithNum() const
{
  G4cout << "Command directory path : " << pathName << G4endl;
  if( guidance != NULL ) guidance->List();
  G4int i = 0;
  G4cout << " Sub-directories : " << G4endl;
  G4int n_treeEntry = tree.size();
  for( G4int i_thTree = 0; i_thTree < n_treeEntry; i_thTree++ )
  {
    i++;
    G4cout << " " << i << ") " << tree[i_thTree]->GetPathName() 
	 << "   " << tree[i_thTree]->GetTitle() << G4endl;
  }
  G4cout << " Commands : " << G4endl;
  G4int n_commandEntry = command.size();
  for( G4int i_thCommand = 0; i_thCommand < n_commandEntry; i_thCommand++ )
  {
    i++;
    G4cout << " " << i << ") " << command[i_thCommand]->GetCommandName() 
	 << " * " << command[i_thCommand]->GetTitle() << G4endl;
  }
}

void G4UIcommandTree::List() const
{
  ListCurrent();
  G4int n_commandEntry = command.size();
  for( G4int i_thCommand = 0; i_thCommand < n_commandEntry; i_thCommand++ )
  {
    command[i_thCommand]->List();
  }
  G4int n_treeEntry = tree.size();
  for( G4int i_thTree = 0; i_thTree < n_treeEntry; i_thTree++ )
  {
    tree[i_thTree]->List();
  }
}

G4String G4UIcommandTree::CreateFileName(const char* pName)
{
  G4String fn = pName;
  G4int idxs;
  while((idxs=fn.index("/"))!=G4int(std::string::npos))
  { fn(idxs) = '_'; }
  fn += ".html";
  return fn;
}

G4String G4UIcommandTree::ModStr(const char* strS)
{
  G4String sx;
  G4String str = strS;
  for(G4int i=0;i<G4int(str.length());i++)
  {
    char c = str(i);
    switch(c)
    {
    case '<':
      sx += "&lt;"; break;
    case '>':
      sx += "&gt;"; break;
    case '&':
      sx += "&amp;"; break;
    default:
      sx += c;
    }
  }
  return sx;
}

void G4UIcommandTree::CreateHTML()
{
  G4String ofileName = CreateFileName(pathName);
  std::ofstream oF(ofileName, std::ios::out);

  oF << "<html><head><title>Commands in " << ModStr(pathName) << "</title></head>" << G4endl;
  oF << "<body bgcolor=\"#ffffff\"><h2>" << ModStr(pathName) << "</h2><p>" << G4endl;

  if( guidance != NULL ) 
  {
    for(G4int i=0;i<guidance->GetGuidanceEntries();i++)
    { oF << ModStr(guidance->GetGuidanceLine(i)) << "<br>" << G4endl; }
  }

  oF << "<p><hr><p>" << G4endl;
  
  oF << "<h2>Sub-directories : </h2><dl>" << G4endl;
  for( G4int i_thTree = 0; i_thTree < G4int(tree.size()); i_thTree++ )
  {
    oF << "<p><br><p><dt><a href=\"" << CreateFileName(tree[i_thTree]->GetPathName())
       << "\">" << ModStr(tree[i_thTree]->GetPathName()) << "</a>" << G4endl;
    oF << "<p><dd>" << ModStr(tree[i_thTree]->GetTitle()) << G4endl;
    tree[i_thTree]->CreateHTML();
  }

  oF << "</dl><p><hr><p>" << G4endl;
  
  oF << "<h2>Commands : </h2><dl>" << G4endl;
  for( G4int i_thCommand = 0; i_thCommand < G4int(command.size()); i_thCommand++ )
  {
    G4UIcommand* cmd = command[i_thCommand];
    oF << "<p><br><p><dt><b>" << ModStr(cmd->GetCommandName());
    if(cmd->GetParameterEntries()>0)
    {
      for(G4int i_thParam=0;i_thParam<cmd->GetParameterEntries();i_thParam++)
      { oF << " [<i>" << ModStr(cmd->GetParameter(i_thParam)->GetParameterName()) << "</i>]"; }
    }
    oF << "</b>" << G4endl;
    oF << "<p><dd>" << G4endl;
    for(G4int i=0;i<cmd->GetGuidanceEntries();i++)
    { oF << ModStr(cmd->GetGuidanceLine(i)) << "<br>" << G4endl; }
    if(!(cmd->GetRange()).isNull())
    { oF << "<p><dd>Range : " << ModStr(cmd->GetRange()) << G4endl; }
    std::vector<G4ApplicationState>* availabelStateList = cmd->GetStateList();
    if(availabelStateList->size()==6)
    { oF << "<p><dd>Available at all Geant4 states." << G4endl; }
    else
    {
      oF << "<p><dd>Available Geant4 state(s) : ";
      for(G4int ias=0;ias<G4int(availabelStateList->size());ias++)
      { oF << G4StateManager::GetStateManager()->GetStateString((*availabelStateList)[ias]) << " " << G4endl; }
    }
    if(cmd->GetParameterEntries()>0)
    {
      oF << "<p><dd>Parameters<table border=1>" << G4endl;
      for(G4int i_thParam=0;i_thParam<cmd->GetParameterEntries();i_thParam++)
      {
        G4UIparameter* prm = cmd->GetParameter(i_thParam);
        oF << "<tr><td>" << ModStr(prm->GetParameterName()) << G4endl;
        oF << "<td>type " << prm->GetParameterType() << G4endl;
        oF << "<td>";
        if(prm->IsOmittable())
        { 
          oF << "Omittable : ";
          if(prm->GetCurrentAsDefault())
          { oF << "current value is used as the default value." << G4endl; }
          else
          { oF << "default value = " << prm->GetDefaultValue() << G4endl; }
        }
        oF << "<td>";
        if(!(prm->GetParameterRange()).isNull())
        { oF << "Parameter range : " << ModStr(prm->GetParameterRange()) << G4endl; }
        else if(!(prm->GetParameterCandidates()).isNull())
        { oF << "Parameter candidates : " << ModStr(prm->GetParameterCandidates()) << G4endl; }
      }
      oF << "</table>" << G4endl;
    }

  }
  
  oF << "</dl></body></html>" << G4endl;
  oF.close();
}


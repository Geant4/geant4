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
//
//
// $Id: G4UIcommandTree.cc,v 1.11 2002-04-26 22:03:35 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4UIcommandTree.hh"
#include "G4StateManager.hh"
#include "g4std/fstream"
#include "G4ios.hh"

G4UIcommandTree::G4UIcommandTree()
:guidance(NULL)
{ }

G4UIcommandTree::G4UIcommandTree(const char * thePathName)
:guidance(NULL)
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

void G4UIcommandTree::AddNewCommand(G4UIcommand *newCommand)
{
  G4String commandPath = newCommand->GetCommandPath();
  G4String remainingPath = commandPath;
  remainingPath.remove(0,pathName.length());
  if( remainingPath.isNull() )
  {
    guidance = newCommand;
    return;
  }
  G4int i = remainingPath.first('/');
  if( i == G4int(G4std::string::npos) )
  {
    // Find command
    G4int n_commandEntry = command.size();
    for( G4int i_thCommand = 0; i_thCommand < n_commandEntry; i_thCommand++ )
    {
      if( remainingPath == command[i_thCommand]->GetCommandName() )
      { return; }
    }
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
	tree[i_thTree]->AddNewCommand( newCommand );
	return; 
      }
    }
    G4UIcommandTree * newTree = new G4UIcommandTree( nextPath );
    tree.push_back( newTree );
    newTree->AddNewCommand( newCommand );
    return;
  }
}

void G4UIcommandTree::RemoveCommand(G4UIcommand *aCommand)
{
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
    if( i == G4int(G4std::string::npos) )
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
    	  if(n_commandRemain==0)
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

G4UIcommand * G4UIcommandTree::FindPath(const char* commandPath)
{
  G4String remainingPath = commandPath;
  if( remainingPath.index( pathName ) == G4std::string::npos )
  { return NULL; }
  remainingPath.remove(0,pathName.length());
  G4int i = remainingPath.first('/');
  if( i == G4int(G4std::string::npos) )
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

void G4UIcommandTree::ListCurrent()
{
  G4cout << "Command directory path : " << pathName << G4endl;
  if( guidance != NULL ) guidance->List();
  G4cout << " Sub-directories : " << G4endl;
  G4int n_treeEntry = tree.size();
  for( G4int i_thTree = 0; i_thTree < n_treeEntry; i_thTree++ )
  {
    G4cout << "   " << tree[i_thTree]->GetPathName() 
	 << "   " << tree[i_thTree]->GetTitle() << G4endl;
  }
  G4cout << " Commands : " << G4endl;
  G4int n_commandEntry = command.size();
  for( G4int i_thCommand = 0; i_thCommand < n_commandEntry; i_thCommand++ )
  {
    G4cout << "   " << command[i_thCommand]->GetCommandName() 
	 << " * " << command[i_thCommand]->GetTitle() << G4endl;
  }
}

void G4UIcommandTree::ListCurrentWithNum()
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

void G4UIcommandTree::List()
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
  while((idxs=fn.index("/"))!=G4int(G4std::string::npos))
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
  G4std::ofstream oF(ofileName, G4std::ios::out);

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
    G4std::vector<G4ApplicationState>* availabelStateList = cmd->GetStateList();
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


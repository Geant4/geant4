// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcommandTree.cc,v 1.5 2001-02-08 06:07:20 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4UIcommandTree.hh"
#include "G4ios.hh"

G4UIcommandTree::G4UIcommandTree()
:guidance(NULL)
{ }

G4UIcommandTree::G4UIcommandTree(G4String thePathName)
:guidance(NULL)
{
  pathName = thePathName;
}

G4UIcommandTree::G4UIcommandTree(const char * thePathName)
:guidance(NULL)
{
  pathName = thePathName;
}

G4UIcommandTree::~G4UIcommandTree()
{
  int i;
  int n_treeEntry = tree.size();
  for( i=0; i < n_treeEntry; i++ )
  { delete tree[i]; }
}

int G4UIcommandTree::operator==(const G4UIcommandTree &right) const
{
  return ( pathName == right.GetPathName() );
}

int G4UIcommandTree::operator!=(const G4UIcommandTree &right) const
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
  int i = remainingPath.first('/');
  if( i == G4std::string::npos )
  {
    // Find command
    int n_commandEntry = command.size();
    for( int i_thCommand = 0; i_thCommand < n_commandEntry; i_thCommand++ )
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
    int n_treeEntry = tree.size();
    for( int i_thTree = 0; i_thTree < n_treeEntry; i_thTree++ )
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
    int i = remainingPath.first('/');
    if( i == G4std::string::npos )
    {
      // Find command
      int n_commandEntry = command.size();
      for( int i_thCommand = 0; i_thCommand < n_commandEntry; i_thCommand++ )
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
      int n_treeEntry = tree.size();
      for( int i_thTree = 0; i_thTree < n_treeEntry; i_thTree++ )
      {
        if( nextPath == tree[i_thTree]->GetPathName() )
        { 
    	  tree[i_thTree]->RemoveCommand( aCommand );
    	  int n_commandRemain = tree[i_thTree]->GetCommandEntry();
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

G4UIcommand * G4UIcommandTree::FindPath(G4String commandPath)
{
  if( commandPath.index( pathName ) == G4std::string::npos )
  { return NULL; }
  G4String remainingPath = commandPath;
  remainingPath.remove(0,pathName.length());
  int i = remainingPath.first('/');
  if( i == G4std::string::npos )
  {
    // Find command
    int n_commandEntry = command.size();
    for( int i_thCommand = 0; i_thCommand < n_commandEntry; i_thCommand++ )
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
    int n_treeEntry = tree.size();
    for( int i_thTree = 0; i_thTree < n_treeEntry; i_thTree++ )
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
  int i = 0;
  G4cout << " Sub-directories : " << G4endl;
  int n_treeEntry = tree.size();
  for( int i_thTree = 0; i_thTree < n_treeEntry; i_thTree++ )
  {
    G4cout << "   " << tree[i_thTree]->GetPathName() 
	 << "   " << tree[i_thTree]->GetTitle() << G4endl;
  }
  G4cout << " Commands : " << G4endl;
  int n_commandEntry = command.size();
  for( int i_thCommand = 0; i_thCommand < n_commandEntry; i_thCommand++ )
  {
    G4cout << "   " << command[i_thCommand]->GetCommandName() 
	 << " * " << command[i_thCommand]->GetTitle() << G4endl;
  }
}

void G4UIcommandTree::ListCurrentWithNum()
{
  G4cout << "Command directory path : " << pathName << G4endl;
  if( guidance != NULL ) guidance->List();
  int i = 0;
  G4cout << " Sub-directories : " << G4endl;
  int n_treeEntry = tree.size();
  for( int i_thTree = 0; i_thTree < n_treeEntry; i_thTree++ )
  {
    i++;
    G4cout << " " << i << ") " << tree[i_thTree]->GetPathName() 
	 << "   " << tree[i_thTree]->GetTitle() << G4endl;
  }
  G4cout << " Commands : " << G4endl;
  int n_commandEntry = command.size();
  for( int i_thCommand = 0; i_thCommand < n_commandEntry; i_thCommand++ )
  {
    i++;
    G4cout << " " << i << ") " << command[i_thCommand]->GetCommandName() 
	 << " * " << command[i_thCommand]->GetTitle() << G4endl;
  }
}

void G4UIcommandTree::List()
{
  ListCurrent();
  int n_commandEntry = command.size();
  for( int i_thCommand = 0; i_thCommand < n_commandEntry; i_thCommand++ )
  {
    command[i_thCommand]->List();
  }
  int n_treeEntry = tree.size();
  for( int i_thTree = 0; i_thTree < n_treeEntry; i_thTree++ )
  {
    tree[i_thTree]->List();
  }
}



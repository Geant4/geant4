// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VBasicShell.cc,v 1.1 1999-01-07 16:09:37 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4VBasicShell.hh"
#include "G4StateManager.hh"
#include "G4UIcommandTree.hh"
#include "G4UIcommand.hh"
#include "G4UIcommandStatus.hh"
#include "G4UImanager.hh"
#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif

G4VBasicShell::G4VBasicShell()
{
  currentDirectory = "/";
}

G4VBasicShell::~G4VBasicShell() 
{;}

G4String G4VBasicShell::ModifyToFullPathCommand(const char* aCommandLine)
{
  G4String rawCommandLine = aCommandLine;
  if(rawCommandLine.isNull()||rawCommandLine(0)=='\0') return rawCommandLine;
  G4String commandLine = rawCommandLine.strip(G4String::both);
  G4String commandString;
  G4String parameterString;
  int i = commandLine.index(" ");
  if( i != RW_NPOS )
  {
    commandString = commandLine(0,i);
    parameterString = " ";
    parameterString += commandLine(i+1,commandLine.length()-(i+1));
  }
  else
  { commandString = commandLine; }

  G4String fullPathCommandLine
    = ModifyPath( commandString )+parameterString;
  return fullPathCommandLine;
}

G4String G4VBasicShell::GetCurrentWorkingDirectory()
{
  return currentDirectory;
}

G4bool G4VBasicShell::ChangeDirectory(const char* newDir)
{
  G4String aNewPrefix = newDir;
  G4String newPrefix = aNewPrefix.strip(G4String::both);
  G4String newDirectory = ModifyPath( newPrefix );
  if( newDirectory( newDirectory.length() - 1 ) != '/' )
  { newDirectory += "/"; }
  if( FindDirectory( newDirectory ) == NULL )
  { return false; }
  currentDirectory = newDirectory;
  return true;
}

G4UIcommandTree* G4VBasicShell::FindDirectory(const char* dirName)
{
  G4String aDirName = dirName;
  G4String theDir = aDirName.strip(G4String::both);
  G4String targetDir = ModifyPath( theDir );
  if( targetDir( targetDir.length()-1 ) != '/' )
  { targetDir += "/"; }
  G4UIcommandTree* comTree = G4UImanager::GetUIpointer()->GetTree();
  if( targetDir == "/" )
  { return comTree; }
  int idx = 1;
  while( idx < targetDir.length()-1 )
  {
    int i = targetDir.index("/",idx);
    comTree = comTree->GetTree(targetDir(0,i+1));
    if( comTree == NULL ) 
    { return NULL; }
    idx = i+1;
  }
  return comTree;
}

G4UIcommand* G4VBasicShell::FindCommand(const char* commandName)
{
  G4String rawCommandLine = commandName;
  G4String commandLine = rawCommandLine.strip(G4String::both);
  G4String commandString;
  int i = commandLine.index(" ");
  if( i != RW_NPOS )
  { commandString = commandLine(0,i); }
  else
  { commandString = commandLine; }

  G4String targetCom = ModifyPath(commandString);
  return G4UImanager::GetUIpointer()->GetTree()->FindPath(targetCom);
}

G4String G4VBasicShell::ModifyPath(G4String tempPath)
{
  G4String newPath = currentDirectory;
  if( tempPath(0) == '/' )   // full path is given
  { newPath = tempPath; }
  else if( tempPath(0) != '.' ) // add current prefix
  { newPath += tempPath; }
  else if( tempPath(0,2) == "./" ) // add current prefix
  { newPath += tempPath(2,tempPath.length()-2); }
  else                       // swim up with ".."
  {
    while( 1 )
    {
      if( tempPath(0,2) == ".." )
      {
        if( newPath != "/" )
        { 
	  G4String tmpString = newPath(0,newPath.length()-1);
          newPath = newPath(0,tmpString.last('/')+1); 
        }
        if( tempPath == ".." || tempPath == "../" )
        { break; }
        tempPath = tempPath(3,tempPath.length()-3);
      }
      else
      {
        newPath += tempPath;
        break;
      }
    }
  }
  return newPath;
}



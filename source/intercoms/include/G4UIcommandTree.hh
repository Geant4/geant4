// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcommandTree.hh,v 1.6 2001-02-08 06:07:18 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4UIcommandTree_h
#define G4UIcommandTree_h 1


#include "G4UIcommand.hh"
#include "globals.hh"
//#include "g4rw/tpordvec.h"
#include "g4std/vector"

// class description:
//
//  This class is exclusively used by G4UImanager for handling the
// tree structure of the commands. The user MUST NOT construct/use
// this class object.

class G4UIcommandTree 
{
  public:
      G4UIcommandTree();
      G4UIcommandTree(G4String thePathName);
      G4UIcommandTree(const char * thePathName);
      ~G4UIcommandTree();
      int operator==(const G4UIcommandTree &right) const;
      int operator!=(const G4UIcommandTree &right) const;

  public:
      void AddNewCommand(G4UIcommand * newCommand);
      void RemoveCommand(G4UIcommand * aCommand);
      G4UIcommand * FindPath(G4String commandPath);
      void List();
      void ListCurrent();
      void ListCurrentWithNum();

  private:
      G4std::vector<G4UIcommand*> command;
      G4std::vector<G4UIcommandTree*> tree;
      G4UIcommand *guidance;
      G4String pathName;

  public:
      inline const G4UIcommand * GetGuidance() const
      { return guidance; };
      inline const G4String GetPathName() const
      { return pathName; };
      inline G4int GetTreeEntry() const
      { return tree.size(); };
      inline G4int GetCommandEntry() const
      { return command.size(); };
      inline G4UIcommandTree * GetTree(int i)
      { return tree[i-1]; };
      inline G4UIcommandTree * GetTree(G4String comName)
      { 
        for( int i=0; i < tree.size(); i++)
        {
          if( comName == tree[i]->GetPathName() )
          { return tree[i]; }
        }
        return NULL;
      };
      inline G4UIcommand * GetCommand(int i)
      { return command[i-1]; };
      inline const G4String GetTitle() const
      { 
	    if(guidance==NULL)
	    { return G4String("...Title not available..."); }
    	else
	    { return guidance->GetTitle(); }
      };

};

#endif


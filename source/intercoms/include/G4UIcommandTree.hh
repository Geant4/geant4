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
// $Id: G4UIcommandTree.hh 99838 2016-10-07 10:06:37Z gcosmo $
//

#ifndef G4UIcommandTree_h
#define G4UIcommandTree_h 1


#include "G4UIcommand.hh"
#include "globals.hh"
#include <vector>

// class description:
//
//  This class is exclusively used by G4UImanager for handling the
// tree structure of the commands. The user MUST NOT construct/use
// this class object.

class G4UIcommandTree 
{
  public:
      G4UIcommandTree();
      G4UIcommandTree(const char * thePathName);
      ~G4UIcommandTree();
      G4int operator==(const G4UIcommandTree &right) const;
      G4int operator!=(const G4UIcommandTree &right) const;

  public:
      void AddNewCommand(G4UIcommand * newCommand, G4bool workerThreadOnly=false);
      void RemoveCommand(G4UIcommand * aCommand, G4bool workerThreadOnly=false);
      G4UIcommand* FindPath(const char* commandPath) const;
      G4UIcommandTree* FindCommandTree(const char* commandPath);
      G4String CompleteCommandPath(const G4String& commandPath);
      G4String GetFirstMatchedString(const G4String&,const G4String&) const;
      // Complete most available caracters in common into command path in the command line
      // given

      void List() const;
      void ListCurrent() const;
      void ListCurrentWithNum() const;
      void CreateHTML();

  private:
      G4String CreateFileName(const char* pName);
      G4String ModStr(const char* strS);

      std::vector<G4UIcommand*> command;
      std::vector<G4UIcommandTree*> tree;
      G4UIcommand *guidance;
      G4String pathName;
      G4bool broadcastCommands;

  public:
      inline const G4UIcommand * GetGuidance() const
      { return guidance; };
      inline const G4String GetPathName() const
      { return pathName; };
      inline G4int GetTreeEntry() const
      { return tree.size(); };
      inline G4int GetCommandEntry() const
      { return command.size(); };
      inline G4UIcommandTree * GetTree(G4int i)
      { return tree[i-1]; };
      inline G4UIcommandTree * GetTree(const char* comNameC)
      { 
        G4String comName = comNameC;
        for( size_t i=0; i < tree.size(); i++)
        {
          if( comName == tree[i]->GetPathName() )
          { return tree[i]; }
        }
        return NULL;
      };
      inline G4UIcommand * GetCommand(G4int i)
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


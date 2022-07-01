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
// G4UIcommandTree
//
// Class description:
//
// This class is exclusively used by G4UImanager for handling the
// tree structure of the commands. The user MUST NOT construct/use
// this class object

// Author: Makoto Asai (SLAC), 1998
// --------------------------------------------------------------------
#ifndef G4UIcommandTree_hh
#define G4UIcommandTree_hh 1

#include <vector>

#include "globals.hh"
#include "G4UIcommand.hh"

class G4UIcommandTree
{
  public:

    G4UIcommandTree() = default;
    G4UIcommandTree(const char* thePathName);

    ~G4UIcommandTree();

    G4bool operator==(const G4UIcommandTree& right) const;
    G4bool operator!=(const G4UIcommandTree& right) const;

    void AddNewCommand(G4UIcommand* newCommand, G4bool workerThreadOnly = false);
    void RemoveCommand(G4UIcommand* aCommand, G4bool workerThreadOnly = false);
    G4UIcommand* FindPath(const char* commandPath) const;
    G4UIcommandTree* FindCommandTree(const char* commandPath);
    G4String GetFirstMatchedString(const G4String&, const G4String&) const;

    G4String CompleteCommandPath(const G4String& commandPath);
      // Complete most available characters in common into command path in the
      // command line given

    void List() const;
    void ListCurrent() const;
    void ListCurrentWithNum() const;
    void CreateHTML(const G4String& = "");

    inline const G4UIcommand* GetGuidance() const { return guidance; }
    inline const G4String& GetPathName() const { return pathName; }
    inline G4int GetTreeEntry() const { return G4int(tree.size()); }
    inline G4int GetCommandEntry() const { return G4int(command.size()); }
    inline G4UIcommandTree* GetTree(G4int i) { return tree[i - 1]; }
    G4UIcommandTree* GetTree(const char* comNameC);
    inline G4UIcommand* GetCommand(G4int i) { return command[i - 1]; }
    inline const G4String GetTitle() const
    {
      return (guidance == nullptr) ? G4String("...Title not available...")
                                   : guidance->GetTitle();
    }

  private:

    G4String CreateFileName(const char* pName);
    G4String ModStr(const char* strS);

    std::vector<G4UIcommand*> command;
    std::vector<G4UIcommandTree*> tree;
    G4UIcommand* guidance = nullptr;
    G4String pathName;
    G4bool broadcastCommands = true;
    G4bool ifSort = false;
    G4int createHTMLTreeLevel = 0;
};

#endif

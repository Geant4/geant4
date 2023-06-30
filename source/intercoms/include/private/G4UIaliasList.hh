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
// G4UIaliasList
//
// Class description:
//
// This class is exclusively used by G4UImanager for handling the
// alias list.

// Author: M.Asai, 1 October 2001
// --------------------------------------------------------------------
#ifndef G4UIaliasList_hh
#define G4UIaliasList_hh 1

#include "globals.hh"

#include <map>

class G4UIaliasList
{
  public:
    // Add or update alias from aliasName to command path aliasValue
    void ChangeAlias(const char* aliasName, const char* aliasValue);

    // Remove aliasName from list of aliases.
    void RemoveAlias(const char* aliasName);

    // Return command corresponding to alias, or nullptr if alias does not exist.
    // Ownership of the pointer is retained by G4UIaliasList.
    const G4String* FindAlias(const char* aliasName) const;

    // Print alias : command pairs to G4cout
    void List() const;

  private:
    // Adds a new alias/value pair to the map
    void AddNewAlias(const char* aliasName, const char* aliasValue);

    std::map<G4String, G4String> aliases;
};

#endif

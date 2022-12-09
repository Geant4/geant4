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
// Author: M.Asai, 1 October 2001
// --------------------------------------------------------------------

#include "G4UIaliasList.hh"
#include "G4ios.hh"

// --------------------------------------------------------------------
G4UIaliasList::~G4UIaliasList()
{
  std::size_t n_treeEntry = alias.size();
  for(std::size_t i = 0; i < n_treeEntry; ++i)
  {
    delete alias[i];
    delete value[i];
  }
}

// --------------------------------------------------------------------
G4bool G4UIaliasList::operator==(const G4UIaliasList& right) const
{
  return (this == &right);
}

// --------------------------------------------------------------------
G4bool G4UIaliasList::operator!=(const G4UIaliasList& right) const
{
  return (this != &right);
}

// --------------------------------------------------------------------
void G4UIaliasList::AddNewAlias(const char* aliasName, const char* aliasValue)
{
  if(FindAlias(aliasName) != nullptr)
  {
    G4cerr << "Alias <" << aliasName << "> already exists. Command ignored."
           << G4endl;
    return;
  }
  G4String* newAlias = new G4String(aliasName);
  alias.push_back(newAlias);
  G4String* newValue = new G4String(aliasValue);
  value.push_back(newValue);
}

// --------------------------------------------------------------------
void G4UIaliasList::RemoveAlias(const char* aliasName)
{
  G4int i = FindAliasID(aliasName);
  if(i < 0)
  {
    G4cerr << "Alias <" << aliasName << "> does not exist. Command ignored."
           << G4endl;
    return;
  }
  alias.erase(alias.begin() + i);
  value.erase(value.begin() + i);
}

// --------------------------------------------------------------------
void G4UIaliasList::ChangeAlias(const char* aliasName, const char* aliasValue)
{
  G4int i = FindAliasID(aliasName);
  if(i < 0)
  {
    AddNewAlias(aliasName, aliasValue);
    return;
  }
  *(value[i]) = aliasValue;
}

// --------------------------------------------------------------------
G4String* G4UIaliasList::FindAlias(const char* aliasName)
{
  G4int i = FindAliasID(aliasName);
  if(i < 0)
  {
    return nullptr;
  }
  return value[i];
}

// --------------------------------------------------------------------
G4int G4UIaliasList::FindAliasID(const char* aliasName)
{
  std::size_t i_entry = alias.size();
  for(std::size_t i = 0; i < i_entry; ++i)
  {
    if(*(alias[i]) == aliasName)
    {
      return (G4int)i;
    }
  }
  return -1;
}

// --------------------------------------------------------------------
void G4UIaliasList::List()
{
  std::size_t i_entry = alias.size();
  for(std::size_t i1 = 0; i1 < i_entry - 1; ++i1)
  {
    for(std::size_t i2 = i1 + 1; i2 < i_entry; ++i2)
    {
      if(*(alias[i1]) > *(alias[i2]))
      {
        G4String* tmp = alias[i1];
        alias[i1]     = alias[i2];
        alias[i2]     = tmp;
        tmp           = value[i1];
        value[i1]     = value[i2];
        value[i2]     = tmp;
      }
    }
  }

  for(std::size_t i = 0; i < i_entry; ++i)
  {
    G4cout << "  " << *(alias[i]) << " : " << *(value[i]) << G4endl;
  }
}

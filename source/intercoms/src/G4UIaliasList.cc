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
// $Id: G4UIaliasList.cc,v 1.3 2001-10-05 00:50:41 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4UIaliasList.hh"
#include "G4ios.hh"

G4UIaliasList::G4UIaliasList()
{ }

G4UIaliasList::~G4UIaliasList()
{
  int i;
  int n_treeEntry = alias.size();
  for( i=0; i < n_treeEntry; i++ )
  { delete alias[i]; 
    delete value[i]; }
}

int G4UIaliasList::operator==(const G4UIaliasList &right) const
{
  return ( this == &right );
}

int G4UIaliasList::operator!=(const G4UIaliasList &right) const
{
  return ( this != &right );
}

void G4UIaliasList::AddNewAlias(G4String aliasName, G4String aliasValue)
{
  if(FindAlias(aliasName))
  {
    G4cerr << "Alias <" << aliasName << "> already exist. Command ignored."
           << G4endl;
    return;
  }
  G4String* newAlias = new G4String(aliasName);
  alias.push_back(newAlias);
  G4String* newValue = new G4String(aliasValue);
  value.push_back(newValue);
}

void G4UIaliasList::RemoveAlias(G4String aliasName)
{
  G4int i = FindAliasID(aliasName);
  if(i<0)
  {
    G4cerr << "Alias <" << aliasName << "> does not exist. Command ignored."
           << G4endl;
    return;
  }
  alias.erase(alias.begin()+i);
  value.erase(value.begin()+i);
}

void G4UIaliasList::ChangeAlias(G4String aliasName, G4String aliasValue)
{
  G4int i = FindAliasID(aliasName);
  if(i<0)
  {
    AddNewAlias(aliasName,aliasValue);
    return;
  }
  *(value[i]) = aliasValue;
}

G4String* G4UIaliasList::FindAlias(G4String aliasName)
{
  G4int i = FindAliasID(aliasName);
  if(i<0)
  { return 0; }
  return value[i];
}

G4int G4UIaliasList::FindAliasID(G4String aliasName)
{
  G4int i_entry = alias.size();
  for(G4int i=0;i<i_entry;i++)
  { if(*(alias[i])==aliasName) return i; }
  return -1;
}

void G4UIaliasList::List() 
{
  G4int i_entry = alias.size();
  for(G4int i1=0;i1<i_entry-1;i1++)
  for(G4int i2=i1+1;i2<i_entry;i2++)
  {
    if(*(alias[i1])>*(alias[i2]))
    {
      G4String* tmp = alias[i1];
      alias[i1] = alias[i2];
      alias[i2] = tmp;
      tmp = value[i1];
      value[i1] = value[i2];
      value[i2] = tmp;
    }
  }

  for(G4int i=0;i<i_entry;i++)
  { G4cout << "  " << *(alias[i]) << " : " << *(value[i]) << G4endl; }
}


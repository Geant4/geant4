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
// $Id: G4UIaliasList.hh,v 1.4 2002-04-26 22:03:34 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4UIaliasList_h
#define G4UIaliasList_h 1


#include "globals.hh"
#include "g4std/vector"

// class description:
//
//  This class is exclusively used by G4UImanager for handling the
// alias list.
// 

class G4UIaliasList 
{
  public:
      G4UIaliasList();
      ~G4UIaliasList();

  private:
      G4int operator==(const G4UIaliasList &right) const;
      G4int operator!=(const G4UIaliasList &right) const;

  public:
      void RemoveAlias(const char* aliasName);
      void ChangeAlias(const char* aliasName, const char* aliasValue);
      G4String* FindAlias(const char* aliasName);
      void List();

  private:
      void AddNewAlias(const char* aliasName, const char* aliasValue);
      G4int FindAliasID(const char* aliasName);

  private:
      G4std::vector<G4String*> alias;
      G4std::vector<G4String*> value;

};

#endif


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
// $Id: G4UIaliasList.hh,v 1.1 2001-09-30 04:12:51 asaim Exp $
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
      int operator==(const G4UIaliasList &right) const;
      int operator!=(const G4UIaliasList &right) const;

  public:
      void AddNewAlias(G4String aliasName, G4String aliasValue);
      void RemoveAlias(G4String aliasName);
      void ChangeAlias(G4String aliasName, G4String aliasValue);
      G4String* FindAlias(G4String aliasName);
      void List();

  private:
      G4int FindAliasID(G4String aliasName);

  private:
      G4std::vector<G4String*> alias;
      G4std::vector<G4String*> value;

};

#endif


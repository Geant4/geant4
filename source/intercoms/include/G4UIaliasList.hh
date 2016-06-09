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
// $Id: G4UIaliasList.hh,v 1.5 2003/06/16 16:55:30 gunter Exp $
// GEANT4 tag $Name: geant4-05-02 $
//

#ifndef G4UIaliasList_h
#define G4UIaliasList_h 1


#include "globals.hh"
#include <vector>

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
      std::vector<G4String*> alias;
      std::vector<G4String*> value;

};

#endif


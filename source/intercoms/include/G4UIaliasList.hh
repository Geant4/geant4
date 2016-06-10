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
// $Id: G4UIaliasList.hh 67965 2013-03-13 09:35:29Z gcosmo $
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


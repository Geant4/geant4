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


#include "G4DNAReactionTypeManager.hh"
#include "G4VReactionType.hh"

G4DNAReactionTypeManager::G4DNAReactionTypeManager()
    : G4VReactionTypeManager()
{}

G4DNAReactionTypeManager::~G4DNAReactionTypeManager()
{
    Clear();
}

void G4DNAReactionTypeManager::Clear()
{
    for(auto it : fReactionTypeTable)
    {
        if(it.second != nullptr)
        {
            delete it.second;
            it.second = nullptr;
        }
    }
    fReactionTypeTable.clear();
}

void G4DNAReactionTypeManager::SetReactionTypeTable(G4VReactionType* process) 
{
    G4int index = fReactionTypeTable.size();
    fReactionTypeTable[index] = process;
}

const ReactionTypeTable* G4DNAReactionTypeManager::GetReactionTypeTable() const
{
    return &fReactionTypeTable;
}

void G4DNAReactionTypeManager::SetTypeTableByID(std::map<ReactionID, ReactionType> byIDMap)
{
    fReactionTypeByID = byIDMap;
}

ReactionType G4DNAReactionTypeManager::GetReactionTypeByID(ReactionID iD)
{
    return fReactionTypeByID[iD];
}

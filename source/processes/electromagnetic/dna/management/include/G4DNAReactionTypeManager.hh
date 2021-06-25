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

#ifndef G4DNAReactionTypeManager_hh
#define G4DNAReactionTypeManager_hh 1

#include <map>
#include "G4VReactionType.hh"

using ReactionTypeTable = std::map<G4int , G4VReactionType*>;
using ReactionID = G4int;

class G4VReactionTypeManager
{
public:
    G4VReactionTypeManager() = default;
    virtual ~G4VReactionTypeManager() = default;
    virtual void SetReactionTypeTable(G4VReactionType*) = 0;
    virtual const ReactionTypeTable* GetReactionTypeTable() const = 0;
    virtual void SetTypeTableByID(std::map<ReactionID, ReactionType> byIDMap) = 0;
};
class G4DNAReactionTypeManager: public G4VReactionTypeManager
{
public:
    G4DNAReactionTypeManager();
    ~G4DNAReactionTypeManager() override;
    void SetReactionTypeTable(G4VReactionType* process) override;

    const ReactionTypeTable* GetReactionTypeTable() const override;

    void SetTypeTableByID(std::map<ReactionID, ReactionType> byIDMap) override;
    ReactionType GetReactionTypeByID(ReactionID iD);

private:
    ReactionTypeTable fReactionTypeTable;
    std::map<ReactionID, ReactionType> fReactionTypeByID;
protected:
    void Clear();
};
#endif



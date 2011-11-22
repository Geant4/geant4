#ifndef G4VReactionModel_
#define G4VReactionModel_

#include "globals.hh"
#include "AddClone_def.hh"

class G4DNAMolecularReactionTable;
class G4Molecule;
class G4Track;

class G4VDNAReactionModel
{
public :
    G4VDNAReactionModel();
    G4VDNAReactionModel(const G4VDNAReactionModel&);
    virtual ~G4VDNAReactionModel();

    /** This macro is defined in AddClone_def **/
    ParentToClone(G4VDNAReactionModel)

    virtual void Initialise(const G4Molecule*, const G4Track&) {;}
    virtual void InitialiseToPrint(const G4Molecule*) = 0 ;
    virtual G4double GetReactionRadius(const G4Molecule*, const G4Molecule*) = 0;
    virtual G4double GetReactionRadius(const int) = 0;
    virtual G4bool FindReaction(const G4Track&, const G4Track&,
                                const G4double /*reactionRadius*/,
                                G4double& /*separationDistance*/,  // To be calculated
                                const G4bool /*hasReachedUserTimeLimit*/) = 0;

    inline void SetReactionTable(const G4DNAMolecularReactionTable*);
    inline const G4DNAMolecularReactionTable* GetReactionTable();

protected :
    G4VDNAReactionModel& operator=(const G4VDNAReactionModel&);
    const G4DNAMolecularReactionTable* fReactionTable ;
};

inline void G4VDNAReactionModel::SetReactionTable(const G4DNAMolecularReactionTable* table)
{
    fReactionTable = table ;
}

inline const G4DNAMolecularReactionTable* G4VDNAReactionModel::GetReactionTable()
{
    return fReactionTable ;
}
#endif

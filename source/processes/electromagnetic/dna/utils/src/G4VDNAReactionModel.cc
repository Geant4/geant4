#include "G4VDNAReactionModel.hh"
#include "G4DNAMolecularReactionTable.hh"

G4VDNAReactionModel::G4VDNAReactionModel()
{
    fReactionTable = 0 ;
}

G4VDNAReactionModel::G4VDNAReactionModel(const G4VDNAReactionModel& right)
{
    fReactionTable = right.fReactionTable ;
}

G4VDNAReactionModel::~G4VDNAReactionModel()
{
    fReactionTable = 0;
}

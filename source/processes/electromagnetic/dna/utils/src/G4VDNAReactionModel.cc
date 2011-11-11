#include "G4VDNAReactionModel.hh"

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

G4VDNAReactionModel& G4VDNAReactionModel::operator=(const G4VDNAReactionModel& right)
{
    if(this == &right) return *this;
    fReactionTable = right.fReactionTable ;
    return *this;
}

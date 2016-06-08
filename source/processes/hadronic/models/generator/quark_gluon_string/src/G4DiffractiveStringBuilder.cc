#include "G4DiffractiveStringBuilder.hh"

//***************************************************************************************************

G4DiffractiveStringBuilder::G4DiffractiveStringBuilder()
   {
   }

G4DiffractiveStringBuilder::G4DiffractiveStringBuilder(const G4DiffractiveStringBuilder &right)
   {
   }

G4DiffractiveStringBuilder::~G4DiffractiveStringBuilder()
   {
   }

//***************************************************************************************************

G4ExcitedString* G4DiffractiveStringBuilder::BuildString(G4PartonPair * aPair)       
    {
    return  new G4ExcitedString(aPair->GetParton1(), aPair->GetParton2(), aPair->GetDirection());
    }

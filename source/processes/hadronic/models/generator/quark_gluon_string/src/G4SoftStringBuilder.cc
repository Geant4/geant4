#include "G4SoftStringBuilder.hh"

//***************************************************************************************************

G4SoftStringBuilder::G4SoftStringBuilder()
   {
   }

G4SoftStringBuilder::G4SoftStringBuilder(const G4SoftStringBuilder &right)
   {
   }

G4SoftStringBuilder::~G4SoftStringBuilder()
   {
   }

//***************************************************************************************************
       
G4ExcitedString* G4SoftStringBuilder::BuildString(G4PartonPair * aPair)       
    {
    return  new G4ExcitedString(aPair->GetParton1(), aPair->GetParton2(), aPair->GetDirection());
    }
       
//***********************************************************************************************

#include "G4eInverseBremsstrahlung.hh"
#include "G4VEmAdjointModel.hh"
#include "G4AdjointBremsstrahlungModel.hh"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
G4eInverseBremsstrahlung::G4eInverseBremsstrahlung(G4bool whichScatCase,G4String process_name,G4AdjointBremsstrahlungModel* aBremAdjointModel):
				G4VAdjointReverseReaction(process_name,whichScatCase)
{theAdjointEMModel = aBremAdjointModel;
 theAdjointEMModel->SetSecondPartOfSameType(false);
 if (IsScatProjToProjCase) SetIntegralMode(true);
 else   SetIntegralMode(false);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
G4eInverseBremsstrahlung::~G4eInverseBremsstrahlung(){
}

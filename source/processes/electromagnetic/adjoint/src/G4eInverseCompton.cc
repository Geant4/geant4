///////////////////////////////////////////////////////
// File name:     G4eInverseCompton
//
// Author:        Laurent Desorgher
//
// Creation date: 20.11.2006
//
///////////////////////////////////////////////////////
#include "G4eInverseCompton.hh"
#include "G4VEmAdjointModel.hh"
#include "G4AdjointComptonModel.hh"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
G4eInverseCompton::G4eInverseCompton(G4bool whichScatCase,G4String process_name,G4AdjointComptonModel* aComptonAdjointModel):
				G4VAdjointInverseScattering(process_name,whichScatCase)
{theAdjointEMModel = aComptonAdjointModel;
 theAdjointEMModel->SetSecondPartOfSameType(false); 
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
G4eInverseCompton::~G4eInverseCompton(){
}

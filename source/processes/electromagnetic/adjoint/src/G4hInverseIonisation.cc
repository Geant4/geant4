///////////////////////////////////////////////////////
// File name:     G4hInverseIonisation
//
// Author:        Laurent Desorgher
//
// Creation date: 15.02.2009
//
///////////////////////////////////////////////////////
#include "G4hInverseIonisation.hh"
#include "G4VEmAdjointModel.hh"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
G4hInverseIonisation::G4hInverseIonisation(G4bool whichScatCase,G4String process_name,G4AdjointhIonisationModel* aEmAdjointModel):
				G4VAdjointReverseReaction(process_name,whichScatCase)
{theAdjointEMModel = aEmAdjointModel;
 theAdjointEMModel->SetSecondPartOfSameType(false);
 SetIntegralMode(true); 
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
G4hInverseIonisation::~G4hInverseIonisation(){
}

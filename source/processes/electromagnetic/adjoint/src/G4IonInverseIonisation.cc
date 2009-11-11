///////////////////////////////////////////////////////
// File name:     G4IonInverseIonisation
//
// Author:        Laurent Desorgher
//
// Creation date: 25.08.2009
//
///////////////////////////////////////////////////////
#include "G4IonInverseIonisation.hh"
#include "G4VEmAdjointModel.hh"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
G4IonInverseIonisation::G4IonInverseIonisation(G4bool whichScatCase,G4String process_name,G4AdjointIonIonisationModel* aEmAdjointModel):
				G4VAdjointReverseReaction(process_name,whichScatCase)
{theAdjointEMModel = aEmAdjointModel;
 theAdjointEMModel->SetSecondPartOfSameType(false);
 SetIntegralMode(true); 
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
G4IonInverseIonisation::~G4IonInverseIonisation(){
}

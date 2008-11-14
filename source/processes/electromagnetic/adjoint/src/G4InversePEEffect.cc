#include "G4InversePEEffect.hh"
#include "G4VEmAdjointModel.hh"
#include "G4AdjointPhotoElectricModel.hh"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
G4InversePEEffect::G4InversePEEffect(G4String process_name,G4AdjointPhotoElectricModel* aModel):
				G4VAdjointInverseScattering(process_name,false)
{theAdjointEMModel = aModel;
 theAdjointEMModel->SetSecondPartOfSameType(false); 
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
G4InversePEEffect::~G4InversePEEffect(){
}

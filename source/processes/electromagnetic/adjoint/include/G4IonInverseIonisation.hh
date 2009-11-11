/////////////////////////////////////////////////////////////////////////////////
//      Module:		G4IonInverseIonisation
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	25 August 2009 creation by L. Desorgher  		
//
//-------------------------------------------------------------
//	Documentation:
//		Adjoint/reverse discrete ionisation for ions
//

#ifndef G4IonInverseIonisation_h
#define G4IonInverseIonisation_h 1

#include "G4VAdjointReverseReaction.hh"
#include "globals.hh"
#include "G4eIonisation.hh"
#include "G4AdjointIonIonisationModel.hh"
class G4IonInverseIonisation: public G4VAdjointReverseReaction

{
public:

  G4IonInverseIonisation(G4bool whichScatCase, G4String process_name, G4AdjointIonIonisationModel* aEmAdjointModel);
  ~G4IonInverseIonisation();
  
private:
    
};

#endif

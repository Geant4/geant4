/////////////////////////////////////////////////////////////////////////////////
//      Module:		G4eInverseIonisation.hh
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	15 April 2007 creation by L. Desorgher  		
//
//-------------------------------------------------------------
//	Documentation:
//		Adjoint/revrese discrete ionisation
//

#ifndef G4eInverseIonisation_h
#define G4eInverseIonisation_h 1

#include "G4VAdjointReverseReaction.hh"
#include "globals.hh"
#include "G4eIonisation.hh"
#include "G4VEmAdjointModel.hh"
class G4eInverseIonisation: public G4VAdjointReverseReaction

{
public:

  G4eInverseIonisation(G4bool whichScatCase, G4String process_name, G4VEmAdjointModel* aEmAdjointModel);
  ~G4eInverseIonisation();
  
private:
    
};

#endif

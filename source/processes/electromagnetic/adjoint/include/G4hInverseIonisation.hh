/////////////////////////////////////////////////////////////////////////////////
//      Module:		G4hInverseIonisation.hh
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	15 February 2009 creation by L. Desorgher  		
//
//-------------------------------------------------------------
//	Documentation:
//		Adjoint/reverse discrete ionisation for proton
//

#ifndef G4hInverseIonisation_h
#define G4hInverseIonisation_h 1

#include "G4VAdjointReverseReaction.hh"
#include "globals.hh"
#include "G4eIonisation.hh"
#include "G4AdjointhIonisationModel.hh"
class G4hInverseIonisation: public G4VAdjointReverseReaction

{
public:

  G4hInverseIonisation(G4bool whichScatCase, G4String process_name, G4AdjointhIonisationModel* aEmAdjointModel);
  ~G4hInverseIonisation();
  
private:
    
};

#endif

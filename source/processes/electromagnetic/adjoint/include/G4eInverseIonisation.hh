/////////////////////////////////////////////////////////////////////////////////
//      Module:		G4eInverseIonisation.hh
//	Author:       	L. Desorgher
//	Date:		15 April 2007
// 	Organisation: 	SpaceIT GmbH
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

#include "G4VAdjointInverseScattering.hh"
#include "globals.hh"
#include "G4eIonisation.hh"
#include "G4VEmAdjointModel.hh"
class G4eInverseIonisation: public G4VAdjointInverseScattering

{
public:

  G4eInverseIonisation(G4bool whichScatCase, G4String process_name, G4VEmAdjointModel* aEmAdjointModel);
  ~G4eInverseIonisation();
  
private:
    
};

#endif

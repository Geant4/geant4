/////////////////////////////////////////////////////////////////////////////////
//      Module:		G4eInverseCompton.hh
//	Author:       	L. Desorgher
//	Date:		25 October 2007
// 	Organisation: 	SpaceIT GmbH
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	25 October 2007 creation by L. Desorgher  		
//
//-------------------------------------------------------------
//	Documentation:
//		Adjoint/reverse Compton
//


#ifndef G4eInverseCompton_h
#define G4eInverseCompton_h 1

#include "G4VAdjointInverseScattering.hh"
#include "globals.hh"
#include "G4eIonisation.hh"
class G4AdjointComptonModel;
class G4eInverseCompton: public G4VAdjointInverseScattering

{
public:

  G4eInverseCompton(G4bool whichScatCase, G4String process_name, G4AdjointComptonModel* aEmAdjointModel);
  ~G4eInverseCompton();
  
private:
    
};

#endif

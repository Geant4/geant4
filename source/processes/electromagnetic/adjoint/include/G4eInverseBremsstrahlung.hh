/////////////////////////////////////////////////////////////////////////////////
//      Module:		G4eInverseBremstrahlung.hh
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
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
//		Adjoint/reverse bremstrahlung
//


#ifndef G4eInverseBremsstrahlung_h
#define G4eInverseBremsstrahlung_h 1

#include "G4VAdjointReverseReaction.hh"
#include "globals.hh"
#include "G4eIonisation.hh"
class G4AdjointBremsstrahlungModel;
class G4eInverseBremsstrahlung: public G4VAdjointReverseReaction

{
public:

  G4eInverseBremsstrahlung(G4bool whichScatCase, G4String process_name, G4AdjointBremsstrahlungModel* aEmAdjointModel);
  ~G4eInverseBremsstrahlung();
  
private:
    
};

#endif

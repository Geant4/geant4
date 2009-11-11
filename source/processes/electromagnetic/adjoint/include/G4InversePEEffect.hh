/////////////////////////////////////////////////////////////////////////////////
//      Module:		G4InversePEEffect.hh
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
//		Adjoint/reverse photo electric process
//

#ifndef G4InversePEEffect_h
#define G4InversePEEffect_h 1

#include "G4VAdjointReverseReaction.hh"
#include "globals.hh"
class G4AdjointPhotoElectricModel;
class G4InversePEEffect: public G4VAdjointReverseReaction

{
public:

  G4InversePEEffect(G4String process_name, G4AdjointPhotoElectricModel* aModel);
  ~G4InversePEEffect();
  
private:
    
};

#endif

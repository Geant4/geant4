/////////////////////////////////////////////////////////////////////////////////
//      Module:		G4AdjointPEEffect.hh
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
//		Adjoint/reverse photo electric process
//

#ifndef G4InversePEEffect_h
#define G4InversePEEffect_h 1

#include "G4VAdjointInverseScattering.hh"
#include "globals.hh"
class G4AdjointPhotoElectricModel;
class G4InversePEEffect: public G4VAdjointInverseScattering

{
public:

  G4InversePEEffect(G4String process_name, G4AdjointPhotoElectricModel* aModel);
  ~G4InversePEEffect();
  
private:
    
};

#endif

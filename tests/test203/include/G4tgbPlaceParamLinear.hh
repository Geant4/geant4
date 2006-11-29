#ifndef G4tgbPlaceParamLinear_H
#define G4tgbPlaceParamLinear_H 1
#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   HsPlaceParamLinear
Author:      P. Arce
Changes:     01/01: creation  
---------------------------------------------------------------------------*/ 

#include "globals.hh"
#include "G4VPVParameterisation.hh"
#include "G4tgbPlaceParam.hh"
class G4VPhysicalVolume;

class G4tgbPlaceParamLinear : public G4tgbPlaceParam
{ 
 public:
  G4tgbPlaceParamLinear( int nCopies, double offset, double step, EAxis axis); 
  virtual ~G4tgbPlaceParamLinear(){};
  void ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const;

 private:

};

#endif



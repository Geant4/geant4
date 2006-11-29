#ifndef G4tgbPlaceParamCircle_H
#define G4tgbPlaceParamCircle_H 1
#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   HsPlaceParamCircle
Author:      P. Arce
Changes:     01/01: creation  
---------------------------------------------------------------------------*/ 

#include "globals.hh"
#include "G4VPVParameterisation.hh"
#include "G4tgbPlaceParam.hh"
class G4VPhysicalVolume;

class G4tgbPlaceParamCircle : public G4tgbPlaceParam
{ 
 public:
  G4tgbPlaceParamCircle( int nCopies, double offset, double step, EAxis axis, double radius); 
  virtual ~G4tgbPlaceParamCircle(){};
  void ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const;

 private:
  double theRadius;

};

#endif



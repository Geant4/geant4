#ifndef G4tgbPlaceParam_H
#define G4tgbPlaceParam_H 1
#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   HsPlaceParam
Author:      P. Arce
Changes:     01/01: creation  
---------------------------------------------------------------------------*/ 

#include "globals.hh"
#include "geomdefs.hh"
#include "G4VPVParameterisation.hh"
#include "G4tgbPlaceParam.hh"
class G4VPhysicalVolume;

class G4tgbPlaceParam : public G4VPVParameterisation
{ 
 public:
  G4tgbPlaceParam( int nCopies, double offset, double step, EAxis axis ): theNoCopies(nCopies), theOffset(offset), theStep(step), theAxis(axis) { };
  virtual ~G4tgbPlaceParam(){};
  void ComputeTransformation(const G4int ,G4VPhysicalVolume *) const { };


 protected:

  int theNoCopies;
  double theOffset;
  double theStep;
  EAxis theAxis;
};

#endif



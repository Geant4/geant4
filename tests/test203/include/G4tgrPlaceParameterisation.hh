#ifndef G4tgrPlaceParameterisation_h
#define G4tgrPlaceParameterisation_h

#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgrPlaceParameterisation
Author:      P. Arce
Changes:     01/01: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Class to descripe the positioning of a G4tgrVolume inside another G4tgrVolume as a parameterised volume. Several types are possible:
- Change only the position and rotation for each copy
- Change also the dimensions
- Change also the solid type
HarpDD does not contain the parameterisations, so the data is just stored in this class, without any calculation of the positions of each copy 
*/
//!:POS_PARAM "volu_name" copyNo "parent_name" "parametrisation_type" number_copies step offset extra_data(n words)

//----------------------------------------------------------------------------  
#include "geomdefs.hh"
#include "G4tgrPlace.hh"
#include "CLHEP/Vector/ThreeVector.h"

class G4tgrPlaceParameterisation : public G4tgrPlace
{

 public:
  G4tgrPlaceParameterisation(){ };
  ~G4tgrPlaceParameterisation(){ };

  // creates an object passing the parameters
  G4tgrPlaceParameterisation( const vector<G4String>& );
 
  //! access functions
  int GetNumberOfCopies() const {return theNoCopies;}
  double GetStep() const {return theStep;}
  double GetOffset() const {return theOffset;}
  G4String GetParamType() const {return theParamType;}
  EAxis GetAxis() const {return theAxis;}
  vector<double> GetExtraData() const {return theExtraData;}

 private:
  int theNoCopies;
  double theStep;
  double theOffset;
  G4String theParamType;
  EAxis theAxis;
  //! extra data needed to build the parametrisation (radius if 'CIRCLE', dimensions if 'CHANGE_DIMENSIONS', ...)
  vector<double> theExtraData;

};
#endif

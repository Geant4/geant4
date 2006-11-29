#ifndef G4tgrPlaceSimple_h
#define G4tgrPlaceSimple_h

#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgrPlaceSimple
Author:      P. Arce
Changes:     01/01: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Class to descripe a simple positioning of a G4tgrVolume inside another G4tgrVolume */ 
//----------------------------------------------------------------------------  

#include <vector>
#include "CLHEP/Vector/ThreeVector.h"
#include "G4tgrPlace.hh"


class G4tgrPlaceSimple : public G4tgrPlace
{
 public:
  G4tgrPlaceSimple(){ };
  virtual ~G4tgrPlaceSimple(){ };

  G4tgrPlaceSimple( const vector<G4String>& wl );

  //acces private data members
 public:
  virtual G4String GetRotMatName() const {return theRotMatName;}
  virtual CLHEP::Hep3Vector GetPlacement() const {return thePlace;}

 protected:
  //! the position with respect to parent
  CLHEP::Hep3Vector thePlace;

  //! the rotation matrix (by name, as the rotations matrices are not yet created)
  G4String theRotMatName;

};
#endif


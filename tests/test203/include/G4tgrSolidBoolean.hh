#ifndef G4tgrSolidBoolean_h
#define G4tgrSolidBoolean_h

#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgrSolidBoolean    
Author:      P. Arce
Changes:     12/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! A G4tgrSolidBoolean is an G4tgrSolid object with a solid made out of the boolean operation of two solids. 
The type of operation can be:  Union, Substraction, Intersection
*/  

#include <vector>
#include "G4tgrSolid.hh"


//---------------------------------------------------------------------------- 
class G4tgrSolidBoolean : public G4tgrSolid { 

 public:    
  G4tgrSolidBoolean(const vector<G4String>& wl);
  virtual ~G4tgrSolidBoolean(){ };


 public:
  //! public access functions

  inline const G4tgrSolid* GetSolid( int ii) const;
  virtual G4String GetRelativeRotMatName() const {return theRelativeRotMatName;}
  virtual G4ThreeVector GetRelativePlace() const {return theRelativePlace;}

//! private data members
 private:   
  //- ! Solid types (Box, Tube, etc) of the solids that compose the boolean solid

  //! Vectors of parameters. 
  vector< vector<double>* > theSolidParams; 

  //! relative placement and rotation of solid 2 w.r.t. solid 1
  G4String theRelativeRotMatName;
  G4ThreeVector theRelativePlace;

  //! the two G4tgrSolid's that combine to make this one
  vector<const G4tgrSolid*> theSolids;

};

#endif
 

inline const G4tgrSolid* G4tgrSolidBoolean::GetSolid( int ii ) const
{  
  if(ii != 0 && ii != 1) {
    cerr << "!!!EXIITING G4tgrSolidBoolean::GetSolid. only two G4tgrSolids (0,1), and you are asking for " << ii << endl;
    exit(1);
  }
  return theSolids[ii];
}

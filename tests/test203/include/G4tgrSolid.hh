#ifndef G4tgrSolid_h
#define G4tgrSolid_h

#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgrSolid    
Author:      P. Arce

---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Class to manage the solid  
*/  

#include <vector>
#include <map>
#include "globals.hh"
#include "G4ThreeVector.hh"
using namespace std;

//----------------------------------------------------------------------------  
class G4tgrSolid { 

 public:    

  G4tgrSolid(){};
  G4tgrSolid( const vector<G4String>& wl);
  virtual ~G4tgrSolid(){ };

  //! Public methods to get and set private data members  
  const G4String GetName() const {return theName;}
  const G4String GetType() const {return theType;}
  const vector< vector<double>* > GetSolidParams() const {return theSolidParams;} 
  virtual G4String GetRelativeRotMatName() const{ return G4String("");};
  virtual G4ThreeVector GetRelativePlace() const{ return G4ThreeVector(0,0,0);};

private:
  void FillSolidParams( const std::vector<G4String>&  wl );

//! private data members
 protected:   
  //! Name of the solid
  G4String theName;   
  //! Type of the solid (Simple, Boolean) //?? needed ??
  G4String theType;   
  //! Type of the solid (Box, Tube, ..., Boolean_UNION, Boolean_SUBS, Boolean_INTERS)

  // Vectors of parameters. 
  vector< vector<double>* > theSolidParams; 

};

#endif

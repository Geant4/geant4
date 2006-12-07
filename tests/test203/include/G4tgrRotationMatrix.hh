#ifndef G4tgrRotationMatrix_h
#define G4tgrRotationMatrix_h

using namespace std;
#include "globals.hh"
#include <vector>
/*---------------------------------------------------------------------------   
ClassName:   G4tgrRotationMatrix    
Author:      P. Arce
Changes:     15/07/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Transient class of a rotation matrix.
*/  

#include <vector>
enum RotMatInputType{rm3,rm6,rm9};

//----------------------------------------------------------------------------  
class G4tgrRotationMatrix { 

 public:    
  G4tgrRotationMatrix(){ };
  ~G4tgrRotationMatrix(){ };

  //! construct the G4tgrRotationMatrix (fill its data members) interpreting the data in the list of words 'wl' 
  G4tgrRotationMatrix( const vector<G4String>& wl );

  G4String GetName() {return theName;}
  std::vector<double>& GetValues() {return theValues;}

 private:
  G4String theName;
  std::vector<double> theValues;
  //thetaX,phiX,thetaY,phiY,thetaZ,phiZ;

  RotMatInputType theInputType;
};

#endif
 

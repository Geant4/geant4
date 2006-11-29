#ifndef G4tgbRotationMatrix_h
#define G4tgbRotationMatrix_h
#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgbRotationMatrix    
Author:      P. Arce
Changes:     15/07/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Transient class of a rotation matrix.
Build a G4RotationMatrix, of each rotation matrix.
*/  

#include <vector>
#include <string>
#include "G4tgrRotationMatrix.hh"
#include "G4RotationMatrix.hh"

//----------------------------------------------------------------------------  
class G4tgbRotationMatrix { 

 public:    
  G4tgbRotationMatrix(){ };
  ~G4tgbRotationMatrix(){ };

  //! construct the G4tgbRotationMatrix (fill its data members) interpreting the data in the list of words 'wl' 
  G4tgbRotationMatrix( G4tgrRotationMatrix* tgr );

  //! vuild a G4RotationMatrix transforming theValues
  G4RotationMatrix* BuildG4RotMatrix( );
  G4RotationMatrix* BuildG4RotMatrixFrom3( std::vector<double>& values );
  G4RotationMatrix* BuildG4RotMatrixFrom6( std::vector<double>& values );
  G4RotationMatrix* BuildG4RotMatrixFrom9( std::vector<double>& values );

  G4String GetName(){
    return theTgrRM->GetName(); 
  }

 private:
  G4tgrRotationMatrix* theTgrRM;
  G4RotationMatrix* theG4RM;

};

#endif
 

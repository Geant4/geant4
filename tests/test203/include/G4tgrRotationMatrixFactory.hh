#ifndef G4tgrRotationMatrixFactory_h
#define G4tgrRotationMatrixFactory_h
#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgrRotationMatrixFactory 
Author:      P. Arce
Changes:     15/07/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Singleton class to manage the building of transient rotation matrix
 */

#include "G4tgrRotationMatrix.hh"
#include <map>

typedef map< G4String, G4tgrRotationMatrix* > mstgrrotm;

class G4tgrRotationMatrixFactory 
{
 public:
  ~G4tgrRotationMatrixFactory();
  
  //! get only instance (it it does not exists, create it)
  static G4tgrRotationMatrixFactory* GetInstance();
  
  //! Build an G4tgrRotationMatrix and add it to theTgrRotMats
  G4tgrRotationMatrix* AddRotMatrix( const vector<G4String>& wl );

  //! Look for an G4tgrRotationMatrix and if not found return 0
  G4tgrRotationMatrix* FindRotMatrix(const G4String& rotm);

  const mstgrrotm& GetRotMatMap() const {return theTgrRotMats;}
  vector<G4tgrRotationMatrix*> GetRotMatList() const {return theTgrRotMatList;}

  //! dump list of rotation matrices
  void DumpRotmList();

 private:
   G4tgrRotationMatrixFactory(){};
   
   static G4tgrRotationMatrixFactory* theInstance;
   
   vector<G4tgrRotationMatrix*> theTgrRotMatList;
   mstgrrotm theTgrRotMats;

};

#endif

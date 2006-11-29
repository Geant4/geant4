#ifndef G4tgbRotationMatrixMgr_h
#define G4tgbRotationMatrixMgr_h
#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgbRotationMatrixMgr 
Author:      P. Arce
Changes:     15/07/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Singleton class to manage the building of transient rotation matrix, as well as the construction of the corresponding G4RotationMatrix's 
 */

#include <map>
#include "G4tgbRotationMatrix.hh"
#include "G4RotationMatrix.hh"

typedef map< G4String, G4tgbRotationMatrix*, less<G4String> > mstgbrotm;
typedef map< G4String, G4RotationMatrix*, less<G4String> > msg4rotm;


class G4tgbRotationMatrixMgr 
{
 public:
  ~G4tgbRotationMatrixMgr();
  
  //! get only instance (it it does not exists, create it)
  static G4tgbRotationMatrixMgr* GetInstance();
  
  void CopyRotMats();

  //! Look for a G4RotationMatrix and if not found create it from the corresponding G4tgbRotationMatrix
  G4RotationMatrix* FindOrBuildG4RotMatrix(const G4String& name);
  //! Look for a G4RotationMatrix and if not found return 0
  G4RotationMatrix* FindG4RotMatrix(const G4String& name);
  
  //! Look for an G4tgbRotationMatrix and if not found exit
  G4tgbRotationMatrix* FindOrBuildTgbRotMatrix(const G4String& name);
  //! Look for an G4tgbRotationMatrix and if not found return 0
  G4tgbRotationMatrix* FindTgbRotMatrix(const G4String& name);


 public:
  const mstgbrotm GetTgbRotMatList() const{
    return theTgbRotMats;
  }
  const msg4rotm& GetG4RotMatList() const { 
    return theG4RotMats;
  }

 private:
   G4tgbRotationMatrixMgr(){};
   
   static G4tgbRotationMatrixMgr* theInstance;
   
 private:
   mstgbrotm theTgbRotMats;
   msg4rotm theG4RotMats;
 
};

ostream& operator<<(ostream&, const G4RotationMatrix &);
#endif

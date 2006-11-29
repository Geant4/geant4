#ifndef G4tgbSensDetMgr_h
#define G4tgbSensDetMgr_h
#include "globals.hh"

/*---------------------------------------------------------------------------   
ClassName:   G4tgbSensDetMgr    
Author:      P. Arce
Changes:     03/08/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Class to manage the G4 sensitive detectors
It is a singleton, accesed always with GetInstance() */ 

#include <vector>
#include <map.h>
#include <multimap.h>

class G4VSensitiveDetector;

//----------------------------------------------------------------------------  
class G4tgbSensDetMgr 
{ 
 public:    

  G4tgbSensDetMgr(){ };
  ~G4tgbSensDetMgr(){ };

  void initialize();

  G4VSensitiveDetector* findSD( const G4String& name );
  void registerSD( const G4String& LVname, G4VSensitiveDetector* sd );
  //! public access functions
 public:
  //! Get the only instance 
  static G4tgbSensDetMgr* GetInstance();  

private:
  static G4tgbSensDetMgr* theInstance;

};

#endif
 

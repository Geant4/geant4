#ifndef G4tgrParameterMgr_h
#define G4tgrParameterMgr_h

using namespace std;
#include "globals.hh"

/*---------------------------------------------------------------------------   
ClassName:   G4tgrParameterMgr    
Author:      P. Arce
Changes:     26/01/01: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Class to manage the parameters.
It is a singleton, accesed always with GetInstance() */ 

#include <map>
#include <vector>

typedef map< G4String, G4String > mapss;

//----------------------------------------------------------------------------  
class G4tgrParameterMgr 
{ 
 public:    

  G4tgrParameterMgr(){ };
  ~G4tgrParameterMgr(){ };


  /// Get the only instance 
  static G4tgrParameterMgr* GetInstance();  
 
  //! add to theParameterList
  void AddParameter( const vector<G4String>& wl, bool mustBeNew = 0 );

  //! find a Parameter with name 'name'. 
  G4String FindParameter( const G4String& name, bool exists = true );

  //! dump list of parameters
  void DumpList();

//! public data access functions
 public:

 private:
  //! map of Parameter's: G4String is the Parameter name, double is the value
  mapss theParameterList;

  static G4tgrParameterMgr* theInstance;

};

#endif

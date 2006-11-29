// G4tgrParameterMgr
#include "G4tgrParameterMgr.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMaterialFactory.hh"
#include "G4tgrRotationMatrixFactory.hh"
#include "G4tgrFileReader.hh"
#include "G4tgrMessenger.hh"

G4tgrParameterMgr* G4tgrParameterMgr::theInstance = 0;

//-------------------------------------------------------------
G4tgrParameterMgr* G4tgrParameterMgr::GetInstance()
{
  if( !theInstance ) {
    theInstance = new G4tgrParameterMgr;
  }
  return theInstance;
}


//-------------------------------------------------------------
void G4tgrParameterMgr::AddParameter( const vector<G4String>& wl, bool mustBeNew  )
{
  //---------- find first if it exists already
  bool existsAlready = 0;
  FindParameter(wl[1], existsAlready);
  if(existsAlready) {
    if( mustBeNew ) {
      cerr << "!!!! EXITING: Parameter already exists: " << wl[1] << G4endl;
      exit(1);
    } else {
      cerr << "!WARNING: Parameter already exists: " << wl[1] << G4endl;
    }
  }
  
  //---------- Check for miminum number of words read 
  G4tgrUtils::CheckWLsize( wl, 3, WLSIZE_EQ, "Parameter::AddParameter");
  
  //----- convert third argument to double, but then store it as string for later use in CLHEP evaluator 
  //  theParameterList[ wl[1] ] = G4tgrUtils::GetFloat( wl[2] );
  float val = G4tgrUtils::GetFloat( wl[2] );
  theParameterList[ wl[1] ] = G4tgrUtils::ftoa( val );

}


//-------------------------------------------------------------
G4String G4tgrParameterMgr::FindParameter( const G4String& name, bool exists )
{
  G4String par = "";
  
  mapss::iterator sdite = theParameterList.find( name );
  if( sdite == theParameterList.end() ) {
    if( exists ) {
      DumpList();
      cerr << "!!!! EXITING: Parameter not found in list: " << name << G4endl;
      exit(1);
    }
  } else {
    exists = 1;
    par = ((*sdite).second);
  }

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4cout << "parameterMgr: parameter found " << name << " = " << par << G4endl; 
#endif
  return par;
}


//-------------------------------------------------------------
void G4tgrParameterMgr::DumpList()
{
  //---------- Dump number of objects of each class
  G4cout << " @@@@@@@@@@@@@@@@@@ Dumping parameter list " << G4endl;
  mapss::const_iterator cite;
  uint NoPos = 0;
  for( cite = theParameterList.begin();cite != theParameterList.end(); cite++ ) {
    G4cout << (*cite).first << " = " << (*cite).second << G4endl;
  }

}

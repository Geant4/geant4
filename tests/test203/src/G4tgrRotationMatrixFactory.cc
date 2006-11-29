#include "G4tgrRotationMatrixFactory.hh"
#include "G4tgrMessenger.hh"
#include "G4tgrUtils.hh"

G4tgrRotationMatrixFactory * G4tgrRotationMatrixFactory::theInstance = 0;

//-------------------------------------------------------------
G4tgrRotationMatrixFactory* G4tgrRotationMatrixFactory::GetInstance()
{
 if( !theInstance ) {
    theInstance = new G4tgrRotationMatrixFactory;
  }
  return theInstance;
}


//-------------------------------------------------------------
G4tgrRotationMatrixFactory::~G4tgrRotationMatrixFactory()
{
  mstgrrotm::iterator cite;
  for( cite = theTgrRotMats.begin(); cite != theTgrRotMats.end(); cite++) {
    delete (*cite).second;
  }
  theTgrRotMats.clear();

}


//-------------------------------------------------------------
G4tgrRotationMatrix* G4tgrRotationMatrixFactory::AddRotMatrix( const vector<G4String>& wl )
{
  //---------- Check for miminum number of words read 
  if( wl.size() != 5 && wl.size() != 8 && wl.size() != 11 ) {
    G4tgrUtils::DumpVS(wl, "!!!! EXITING:  G4tgrRotationMatrixFactory::AddRotMatrix. Line should have 5, 8 or 11 words ");
    G4Exception("");
  }

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4cout << " add rot matrix " << wl[1] << G4endl;
#endif
  //---------- Look if rotation matrix exists
  if( FindRotMatrix( G4tgrUtils::SubQuotes(wl[1]) ) != 0 ) {
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
      G4tgrUtils::DumpVS( wl, "! WARNING: rotation matrix repeated ", cerr );
    /*    G4tgrUtils::DumpVS( wl, "!!!! EXITING: rotation matrix repeated ", cerr );
	  exit(1); */
#endif
  } 
 
  G4tgrRotationMatrix* rotm = new G4tgrRotationMatrix( wl );
  theTgrRotMats[ rotm->GetName() ] =  rotm;
  theTgrRotMatList.push_back( rotm );
 
  return rotm;
}


//-------------------------------------------------------------
G4tgrRotationMatrix* G4tgrRotationMatrixFactory::FindRotMatrix(const G4String& name)
{
  G4tgrRotationMatrix* rotm = 0;

  mstgrrotm::const_iterator cite = theTgrRotMats.find( name );
  if( cite != theTgrRotMats.end() ) { 
    rotm = (*cite).second;
  } 

  return rotm;

}


//-------------------------------------------------------------
void G4tgrRotationMatrixFactory::DumpRotmList()
{
  G4cout << " @@@@@@@@@@@@@@@@ DUMPING G4tgrRotationMatrix's List " << G4endl;
  mstgrrotm::const_iterator cite;
  for(cite = theTgrRotMats.begin(); cite != theTgrRotMats.end(); cite++) {
    G4cout << " ROTM: " << (*cite).second->GetName() << G4endl;
  }
}

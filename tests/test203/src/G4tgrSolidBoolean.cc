#include "G4tgrSolidBoolean.hh"
#include "G4tgrUtils.hh"
#include "G4tgrVolume.hh"
#include "G4tgrVolumeMgr.hh"
#include "G4tgrMessenger.hh"
#include "G4tgrFileReader.hh"

#include "CLHEP/Units/SystemOfUnits.h"

using namespace CLHEP;


//-------------------------------------------------------------
G4tgrSolidBoolean::G4tgrSolidBoolean( const vector<G4String>& wl ) 
{ 
  //:SOLID/:VOLU VOLU UNION/SUBS/INTERS VOLU1 VOLU2 ROTM POSX POSY POSZ

  if( wl.size() != 9 ) {
    G4tgrUtils::DumpVS(wl, "!!!! EXITING: G4tgrSolidBoolean::G4tgrSolidBoolean Line read with less or more than 9 words ");
    exit(1);
  }

  //---------- set name 
  theName = G4tgrUtils::SubQuotes( wl[1] ); 

  G4tgrVolumeMgr* volmgr = G4tgrVolumeMgr::GetInstance();
  const G4tgrSolid* sol1 = volmgr->FindSolid( G4tgrUtils::SubQuotes( wl[3] ));
  if( !sol1 ) sol1 = volmgr->FindVolume( G4tgrUtils::SubQuotes( wl[3] ), 1)->GetSolid();
  const G4tgrSolid* sol2 = volmgr->FindSolid( G4tgrUtils::SubQuotes( wl[4] ));
  if( !sol2 ) sol2 = volmgr->FindVolume( G4tgrUtils::SubQuotes( wl[4] ), 1)->GetSolid();
  theSolids.push_back( sol1 );
  theSolids.push_back( sol2 );

  /* solid parameters are taken from components
  //--------- set solid parameters (sets of solid parameters of two components)
  //--------- create vector<double> of theSolidParams
  G4cout << " sol1 " << sol1 << G4endl;
  vector< vector<double>* > sp = sol1->GetSolidParams();
  G4cout << " sol1 sp " << sp.size() << " " << sol1->GetSolidParams().size() <<  G4endl;
  theSolidParams.push_back( sp[0] );
  sp = sol2->GetSolidParams();
  theSolidParams.push_back( sp[0] );

  G4cout << " boolsoli sp " << this << " " << GetSolidParams().size() << G4endl;
  */

  //---------- set relative placement and rotation matrix
  theRelativeRotMatName = G4tgrUtils::SubQuotes( wl[5] );
  theRelativePlace = Hep3Vector( G4tgrUtils::GetFloat(wl[6]), G4tgrUtils::GetFloat(wl[7]), G4tgrUtils::GetFloat(wl[8]) );

  //---------- set solid type
  G4String wl2 = wl[2];
  for( size_t ii = 0; ii < wl2.length(); ii++ ){
    wl2[ii] = toupper( wl2[ii] );
  }
  theType = "Boolean_" + wl2;

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) { 
    G4cout << "G4tgrSolidBoolean constructed " << this << " " << GetName() << G4endl;
  }
#endif

  G4tgrVolumeMgr::GetInstance()->RegisterMe( this );

}


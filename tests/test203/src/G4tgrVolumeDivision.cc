#include "G4tgrVolumeDivision.hh"
#include "G4tgrUtils.hh"
#include "G4tgrVolumeMgr.hh"
#include "G4tgrPlace.hh"
#include "G4tgrFileReader.hh"
#include "G4tgrPlaceDivRep.hh"
#include "G4tgrMessenger.hh"

#include "CLHEP/Units/SystemOfUnits.h"

mmss G4tgrVolumeDivision::theSupportedAxis;

//-------------------------------------------------------------
G4tgrVolumeDivision::G4tgrVolumeDivision( const vector<G4String>& wl ) 
{
  //wl: NAME PARENT  MATERIAL AXIS STEP/NDIV OFFSET

  G4tgrUtils::CheckWLsize( wl, 6, WLSIZE_GE, "G4tgrVolumeDivision::G4tgrVolumeDivision");
  G4tgrUtils::CheckWLsize( wl, 8, WLSIZE_LE, "G4tgrVolumeDivision::G4tgrVolumeDivision");

  theType = "VOLDivision";

   //:DIV  NAME PARENT MATERIAL AXIS STEP/NDIV OFFSET

  //---------- set name 
  theName = G4tgrUtils::SubQuotes( wl[1] ); 
  
  //---------- set the pointer to the parent DU
  G4String parentName = G4tgrUtils::SubQuotes(wl[2]);
  G4tgrVolumeMgr::GetInstance()->FindVolume( parentName, 1); // to check if it exits

  //---------- initialize G4tgrPlace
  thePlaceDiv = new G4tgrPlaceDivRep();
  thePlaceDiv->SetParentName( parentName );
  thePlaceDiv->SetType("PlaceDivision");
  thePlaceDiv->SetVolume( this ); 

  //---------- set material name
  theMaterialName = G4tgrUtils::SubQuotes( wl[3] );
 
  //----- set axis of replica
  thePlaceDiv->SetAxis( thePlaceDiv->BuildAxis( G4tgrUtils::SubQuotes(wl[4]) ) ); 

  //------ register parent - child 
  G4tgrVolumeMgr::GetInstance()->RegisterParentChild( parentName, thePlaceDiv );
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4cout << " replica  register parent - child " << G4endl;
#endif

  //---------- set if division is given by number of divisions of by width
  G4bool byStep;
  G4String wl0 = wl[0];
  for( size_t ii = 0; ii < wl0.length(); ii++ ){
    wl0[ii] = toupper( wl0[ii] );
  }

  if( wl0 == ":DIV_NDIV" ) {
    thePlaceDiv->SetDivType( DivByNdiv );
    thePlaceDiv->SetNDiv( G4tgrUtils::GetInt( wl[5] ) );
    if( wl.size() == 7 ) thePlaceDiv->SetOffset( G4tgrUtils::GetFloat( wl[6] )*mm );
  } else if( wl0 == ":DIV_WIDTH" ) {
    thePlaceDiv->SetDivType( DivByWidth );
    thePlaceDiv->SetWidth( G4tgrUtils::GetFloat( wl[5] )*mm );
    if( wl.size() == 7 ) thePlaceDiv->SetOffset( G4tgrUtils::GetFloat( wl[6] )*mm );
  } else if( wl0 == ":DIV_NDIV_WIDTH" ) {
    thePlaceDiv->SetDivType( DivByNdivAndWidth );
    thePlaceDiv->SetNDiv( G4tgrUtils::GetInt( wl[5] ) );
    thePlaceDiv->SetWidth( G4tgrUtils::GetFloat( wl[6] )*mm );
    if( wl.size() == 8 ) thePlaceDiv->SetOffset( G4tgrUtils::GetFloat( wl[7] )*mm );
  } else {
    G4Exception("G4tgrVolumeDivision::G4tgrVolumeDivision division type not supported " + wl[0] );
  }

  theVisibility = 1;
  theRGBColour = new double[3];
  uint ii;
  for(ii=0; ii<3; ii++) theRGBColour[ii] = -1.;
}


#include "G4tgrPlaceSimple.hh"
#include "G4tgrUtils.hh"
#include "G4tgrRotationMatrixFactory.hh"

#include "CLHEP/Units/SystemOfUnits.h"

using namespace CLHEP;

//-------------------------------------------------------------
G4tgrPlaceSimple::G4tgrPlaceSimple( const vector<G4String>& wl )
{
  theType = "PlaceSimple";

  //---------- set the parent name
  theParentName = G4tgrUtils::SubQuotes( wl[3] ); 

  //---------- set the copy number
  theCopyNo = G4tgrUtils::GetInt( wl[2] );

  //---------- set the position with respect to parent
  thePlace = Hep3Vector( G4tgrUtils::GetFloat(wl[5])*mm, G4tgrUtils::GetFloat(wl[6])*mm, G4tgrUtils::GetFloat(wl[7])*mm ); 

  //---------- set the rotation matrix name
  //---------- check if it is going to use the last matrix defined
  if ( wl[4] == ":LAST" ) {
    theRotMatName = ( (*( (G4tgrRotationMatrixFactory::GetInstance()->GetRotMatMap()).rbegin() )).second )->GetName();
    //-    G4cout << " theRotMatName " << theRotMatName << G4endl;
  } else {
    theRotMatName = G4tgrUtils::SubQuotes(wl[4]);
  }

}


#include "G4tgrPlaceParameterisation.hh"
#include "G4tgrUtils.hh"
#include "G4tgrRotationMatrixFactory.hh"

#include "CLHEP/Units/SystemOfUnits.h"

//-------------------------------------------------------------
G4tgrPlaceParameterisation::G4tgrPlaceParameterisation( const vector<G4String>& wl )
{
  theType = "PlaceParam";

  //---------- Check for exact number of words read 
  G4tgrUtils::CheckWLsize( wl, 7, WLSIZE_GE, "G4tgrPlaceParameterisation::ConstructVolume");
  
  //---------- the copy No
  theCopyNo = G4tgrUtils::GetInt( wl[2] )-1;

  //---------- set the parent name
  theParentName = G4tgrUtils::SubQuotes( wl[3] ); 

  //---------- set the number of copies
  theNoCopies = G4tgrUtils::GetInt(wl[5]); 

  //---------- set the step
  theStep = G4tgrUtils::GetFloat(wl[6]); 

  //---------- set the offset
  if( wl.size() > 6 ) { 
    theOffset = G4tgrUtils::GetFloat(wl[7]); 
  } else {
    theOffset = 0.;
  }

  //---------- set the extra data 
  theParamType = G4tgrUtils::SubQuotes( wl[4] );
  if( theParamType == "CIRCLE" ) {
    G4tgrUtils::CheckWLsize( wl, 9, WLSIZE_EQ, "G4tgrPlaceParameterisation::ConstructVolume: CIRCLE");
    theExtraData.push_back( G4tgrUtils::GetFloat(wl[8])*mm );
    theAxis = kPhi;
  }else if( theParamType == "LINEAR_X" || theParamType == "LINEAR_Y" || theParamType == "LINEAR_Z" ) {
    G4tgrUtils::CheckWLsize( wl, 8, WLSIZE_EQ, "G4tgrPlaceParameterisation::ConstructVolume: LINEAR");
    if( theParamType == "LINEAR_X") {
      theAxis = kXAxis;
    } else if( theParamType == "LINEAR_Y") {
      theAxis = kYAxis;
    } else if( theParamType == "LINEAR_Z" ) {
      theAxis = kZAxis;
    }
  } else {
    cerr << "!!!EXITING: G4tgrPlaceParameterisation::ConstructVolume. Parameterisation type not implemented yet: " << theParamType << G4endl
	 << " Please contact Pedro.Arce@cern.ch if you require its implementation " << G4endl;
    exit(1);
  }

}


#include "G4tgrRotationMatrix.hh"
#include "G4tgrRotationMatrixFactory.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMessenger.hh"
#include "CLHEP/Units/SystemOfUnits.h"

//-------------------------------------------------------------
G4tgrRotationMatrix::G4tgrRotationMatrix( const vector<G4String>& wl ) 
{
  theName = G4tgrUtils::SubQuotes( wl[1] );

  switch( wl.size() ){
  case 5:
    theInputType = rm3;
    break;
  case 8:
    theInputType = rm6;
    break;
  case 11:
    theInputType = rm9;
    break;
  default:
    G4Exception("G4tgrRotationMatrix: input line must have 5, 8 or 11 words");
    break;
  }

//-------- Fill matrix values
  uint siz = wl.size() - 2;
  for( uint ii = 0; ii < siz; ii++) {
    theValues.push_back( G4tgrUtils::GetFloat( wl[ii+2] ) );
  }
  if( siz != 9 ) {
    for( uint ii = 0; ii < siz; ii++) {
      theValues[ii] *= deg ;
    }
  }

  /*  thetaX = G4tgrUtils::GetFloat( wl[2] );
  phiX = G4tgrUtils::GetFloat( wl[3] );
  thetaY = G4tgrUtils::GetFloat( wl[4] );
  phiY = G4tgrUtils::GetFloat( wl[5] );
  thetaZ = G4tgrUtils::GetFloat( wl[6] );
  phiZ = G4tgrUtils::GetFloat( wl[7] );
  */
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4cout << " G4tgrRotationMatrix created " << theName << G4endl;
#endif

}



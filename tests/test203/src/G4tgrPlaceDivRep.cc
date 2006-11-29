#include "G4tgrPlaceDivRep.hh"
#include "G4tgrUtils.hh"
#include "G4tgrRotationMatrixFactory.hh"

#include "CLHEP/Units/SystemOfUnits.h"


//-------------------------------------------------------------
G4tgrPlaceDivRep::G4tgrPlaceDivRep()
{
  theOffset = 0.;

}


//-------------------------------------------------------------
G4tgrPlaceDivRep::G4tgrPlaceDivRep( const std::vector<G4String>& wl ) 
{
  theDivType = DivByNdivAndWidth;

  //NAME PARENT AXIS NREP WIDTH OFFSET
  G4tgrUtils::CheckWLsize( wl, 6, WLSIZE_GE, "G4tgrPlaceDivRep::G4tgrPlaceDivRep");
  G4tgrUtils::CheckWLsize( wl, 7, WLSIZE_LE, "G4tgrPlaceDivRep::G4tgrPlaceDivRep");

  theParentName = G4tgrUtils::SubQuotes(wl[2]); 
  theAxis = BuildAxis( G4tgrUtils::SubQuotes(wl[3]) ); 
  theNDiv = G4tgrUtils::GetInt( wl[4] );
  theWidth = G4tgrUtils::GetFloat(wl[5])*mm; //check if it is deg
  if( wl.size() == 7 ) {
    theOffset = G4tgrUtils::GetFloat(wl[6])*mm;  //check if it is deg
  } else {
    theOffset = 0.;
  }

}


//-------------------------------------------------------------
EAxis G4tgrPlaceDivRep::BuildAxis( const G4String& axisName ) 
{
  if( axisName == "X" ) {
    return kXAxis;
  } else if( axisName == "Y" ) {
    return kYAxis;
  } else if( axisName == "Z" ) {
    return kZAxis;
  } else if( axisName == "R" ) {
    return kRho;
  } else if( axisName == "PHI" ) {
    return kPhi;
  } else { 
    cerr << "!!EXITING G4tgrVolumeDivision::GetReplicaAxis. Axis type not found: " << axisName << G4endl
	 << " only valid axis: X, Y, Z, R, PHI " << G4endl;
      exit(1);
  }

}

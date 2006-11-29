#include "G4tgbPlaceParamLinear.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4tgrMessenger.hh"

G4tgbPlaceParamLinear::G4tgbPlaceParamLinear( int nCopies, double offset, double step, EAxis axis) :  G4tgbPlaceParam( nCopies, offset, step, axis )
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4cout << "G4tgbPlaceParamLinear: no copies " << theNoCopies << " = " << nCopies << G4endl
	   << " offset " << theOffset << " = " << offset << G4endl
	   << " step " << theStep << " = " << step << G4endl;
#endif
}

void G4tgbPlaceParamLinear::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume *physVol) const
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) {
    G4cout << "G4tgbPlaceParamLinear::ComputeTransformation" << G4endl;
    G4cout << " no copies " << theNoCopies 
	   << " offset " << theOffset 
	   << " step " << theStep << G4endl;
  }
#endif

  G4ThreeVector origin(0.,0.,0.); 
  double posi = theOffset + copyNo*theStep;
  if( theAxis == kXAxis ) {
    origin.setX( posi ); 
  }else if( theAxis == kYAxis ) {
    origin.setY( posi ); 
  }else if( theAxis == kZAxis ) {
    origin.setZ( posi ); 
  }
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4cout << " G4tgbPlaceParamLinear " << copyNo << " pos " << origin << " rot mat " << " axis " << theAxis << G4endl;
#endif
  //----- set traslation and rotation
  physVol->SetTranslation(origin);
}

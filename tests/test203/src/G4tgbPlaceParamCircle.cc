#include "G4tgbPlaceParamCircle.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4tgrUtils.hh"

G4tgbPlaceParamCircle::G4tgbPlaceParamCircle( int nCopies, double offset, double step, EAxis axis, double radius) :  G4tgbPlaceParam( nCopies, offset, step, axis )
{
  theRadius = radius;
  //-  cout << " G4tgbPlaceParamCircle::G4tgbPlaceParamCircle radiu set " << theRadius << endl;
  //-  cout << " no copies " << theNoCopies << " = " << nCopies << endl
  //-       << " offset " << theOffset << " = " << offset << endl
  //-     << " step " << theStep << " = " << step << endl;
}

void G4tgbPlaceParamCircle::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume *physVol) const
{
  //-  cout << "G4tgbPlaceParamCircle::ComputeTransformation" << endl;
  //- cout << " no copies " << theNoCopies 
  //-    << " offset " << theOffset 
  //-     << " step " << theStep << endl;

  G4ThreeVector origin(1,1,1);
  origin.setTheta(acos(0.));
  double posi = theOffset + copyNo*theStep;
    origin.setPhi( posi );
  origin.setMag( theRadius );
  //----- calculate rotation matrix (so that all volumes point to the centre)
  G4RotationMatrix* rm = new G4RotationMatrix;
  rm->rotateZ( -posi );

  //-  cout << " G4tgbPlaceParamCircle " << copyNo << " pos " << origin << " rot mat " << endl;
  //-  Utils::dumprm(*rm, "rm");
  //----- set traslation and rotation
  physVol->SetTranslation(origin);
  physVol->SetRotation(rm);
}

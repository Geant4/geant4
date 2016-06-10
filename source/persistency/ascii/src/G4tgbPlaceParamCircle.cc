//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4tgbPlaceParamCircle.cc 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgbPlaceParamCircle

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgbPlaceParamCircle.hh"
#include "G4tgrPlaceParameterisation.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMessenger.hh"

// -------------------------------------------------------------------------
G4tgbPlaceParamCircle::~G4tgbPlaceParamCircle()
{
}


// -------------------------------------------------------------------------
G4tgbPlaceParamCircle::
G4tgbPlaceParamCircle( G4tgrPlaceParameterisation* tgrParam )
  : G4tgbPlaceParameterisation(tgrParam)
{
  //---- Get translation and rotation 
  if( tgrParam->GetParamType() == "CIRCLE" )
  {
    CheckNExtraData( tgrParam, 7, WLSIZE_EQ, "G4tgbPlaceParamCircle:");
    theCircleAxis = G4ThreeVector( tgrParam->GetExtraData()[4],
                                   tgrParam->GetExtraData()[5],
                                   tgrParam->GetExtraData()[6] );

    G4ThreeVector zaxis(0.,0.,-1.);
    if( zaxis.cross(theCircleAxis).mag() > 1.E-6 )
    {
      theDirInPlane = zaxis.cross(theCircleAxis);
    }
    else
    { 
      theDirInPlane = theCircleAxis.cross(G4ThreeVector(0.,-1.,0.));
    }
    theAxis = kZAxis;
  }
  else
  {
    CheckNExtraData( tgrParam, 4, WLSIZE_EQ, "G4tgbPlaceParamCircle:");
    if( tgrParam->GetParamType() == "CIRCLE_XY" ) {
      theCircleAxis = G4ThreeVector(0.,0.,1.);
      theDirInPlane = G4ThreeVector(1.,0.,0.);
      theAxis = kZAxis;
    } else if( tgrParam->GetParamType() == "CIRCLE_XZ" ) {
      theCircleAxis = G4ThreeVector(0.,1.,0.);
      theDirInPlane = G4ThreeVector(1.,0.,0.);
      theAxis = kYAxis;
    } else if( tgrParam->GetParamType() == "CIRCLE_YZ" ) {
      theCircleAxis = G4ThreeVector(1.,0.,0.);
      theDirInPlane = G4ThreeVector(0.,1.,0.);
      theAxis = kXAxis;
    }
  }

  if( theCircleAxis.mag() == 0. )
  {
    G4Exception("G4tgbPlaceParamCircle::G4tgbPlaceParamCircle()",
                "InvalidSetup", FatalException, "Circle axis is zero !");
  }
  theCircleAxis /= theCircleAxis.mag();

  theAxis = kZAxis;

  theNCopies = G4int(tgrParam->GetExtraData()[0]);
  theStep = tgrParam->GetExtraData()[1];
  theOffset = tgrParam->GetExtraData()[2];
  theRadius = tgrParam->GetExtraData()[3];

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << " G4tgbPlaceParamCircle::G4tgbPlaceParamCircle():" << G4endl
           << " param type " << tgrParam->GetParamType() << G4endl
           << "   no copies - " << theNCopies << G4endl
           << "   step - " << theStep << G4endl
           << "   offset - " << theOffset << G4endl
           << "   radius - " << theRadius << G4endl
           << "   circle axis - " << theCircleAxis << G4endl
           << "   dir in plane - " << theDirInPlane << G4endl;
  }
#endif
}


// -------------------------------------------------------------------------
void G4tgbPlaceParamCircle::
ComputeTransformation(const G4int copyNo, G4VPhysicalVolume *physVol) const
{ 
  G4double posi = theOffset + copyNo*theStep;
  G4ThreeVector origin = theDirInPlane * theRadius;
  origin.rotate( posi, theCircleAxis );

  //----- Calculate rotation matrix (so that all volumes point to the centre)
  G4RotationMatrix rm;
  rm.rotate( -posi, theCircleAxis );

  //----- Set translation and rotation
  physVol->SetTranslation(origin);
  G4RotationMatrix* pvRm = physVol->GetRotation();
  if( pvRm == 0 )
  {
    pvRm = new G4RotationMatrix;
  }
  *pvRm  = *theRotationMatrix * rm;
  physVol->SetRotation(pvRm);
  physVol->SetCopyNo( copyNo );

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 3 )
  {
    G4cout << " G4tgbPlaceParamCircle::ComputeTransformation():" 
	   << physVol->GetName() << G4endl
           << "   no copies - " << theNCopies  << G4endl
           << "   centre - " << origin << G4endl
           << "   rotation-matrix - " << *pvRm << G4endl;
  }
#endif
}

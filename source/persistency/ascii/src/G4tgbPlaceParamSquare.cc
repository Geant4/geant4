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
// $Id: G4tgbPlaceParamSquare.cc 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgbPlaceParamSquare

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgbPlaceParamSquare.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4tgrMessenger.hh"
#include "G4tgrPlaceParameterisation.hh"

// -------------------------------------------------------------------------
G4tgbPlaceParamSquare::~G4tgbPlaceParamSquare()
{
}


// -------------------------------------------------------------------------
G4tgbPlaceParamSquare::
G4tgbPlaceParamSquare( G4tgrPlaceParameterisation* tgrParam )
  : G4tgbPlaceParameterisation(tgrParam)
{
  //---- Get translation and rotation 
  if( tgrParam->GetParamType() == "SQUARE" )
  {
    CheckNExtraData( tgrParam, 12, WLSIZE_EQ, "G4tgbPlaceParamSquare:");
    theDirection1 = G4ThreeVector( tgrParam->GetExtraData()[6],
                                   tgrParam->GetExtraData()[7],
                                   tgrParam->GetExtraData()[8] );
    theDirection2 = G4ThreeVector( tgrParam->GetExtraData()[9],
                                   tgrParam->GetExtraData()[10],
                                   tgrParam->GetExtraData()[11] );
    theAxis = kZAxis;
  }
  else
  {
    CheckNExtraData( tgrParam, 6, WLSIZE_EQ, "G4tgbPlaceParamSquare:");
    if( tgrParam->GetParamType() == "SQUARE_XY" )
    {
      theDirection1 = G4ThreeVector(1.,0.,0.);
      theDirection2 = G4ThreeVector(0.,1.,0.);
      theAxis = kZAxis;
    }
    else if( tgrParam->GetParamType() == "SQUARE_YZ" )
    {
      theDirection1 = G4ThreeVector(0.,1.,0.);
      theDirection2 = G4ThreeVector(0.,0.,1.);
      theAxis = kXAxis;
    }
    else if( tgrParam->GetParamType() == "SQUARE_XZ" )
    {
      theDirection1 = G4ThreeVector(1.,0.,0.);
      theDirection2 = G4ThreeVector(0.,0.,1.);
      theAxis = kYAxis;
    }
  }

  if( theDirection1.mag() == 0. )
  {
    G4Exception("G4tgbPlaceParamSquare::G4tgbPlaceParamSquare()",
                "InvalidSetup", FatalException, "Direction1 is zero !");
  }
  else
  {
    theDirection1 /= theDirection1.mag();
  }
  if( theDirection2.mag() == 0. )
  {
    G4Exception("G4tgbPlaceParamSquare::G4tgbPlaceParamSquare()",
                "InvalidSetup", FatalException, "Direction2 is zero !");
  }
  else
  {
    theDirection2 /= theDirection2.mag();
  }
  
  theNCopies1 = G4int(tgrParam->GetExtraData()[0]);
  theNCopies2 = G4int(tgrParam->GetExtraData()[1]);
  theStep1 = tgrParam->GetExtraData()[2];
  theStep2 = tgrParam->GetExtraData()[3];
  theOffset1 = tgrParam->GetExtraData()[4];
  theOffset2 = tgrParam->GetExtraData()[5];

  theNCopies = theNCopies1 * theNCopies2;
  theTranslation = theOffset1*theDirection1 + theOffset2*theDirection2;
  
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4cout << "G4tgbPlaceParamSquare: no copies "
           << theNCopies << " = " << theNCopies1
           << " X " << theNCopies2 << G4endl
           << " offset1 " << theOffset1 << G4endl
           << " offset2 " << theOffset1 << G4endl
           << " step1 " << theStep1 << G4endl
           << " step2 " << theStep2 << G4endl
           << " direction1 " << theDirection1 << G4endl
           << " direction2 " << theDirection2 << G4endl
           << " translation " << theTranslation << G4endl;
#endif
}


// -------------------------------------------------------------------------
void G4tgbPlaceParamSquare::
ComputeTransformation(const G4int copyNo, G4VPhysicalVolume *physVol) const
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 3 )
  {
    G4cout << " G4tgbPlaceParamSquare::ComputeTransformation():" 
	   << physVol->GetName() << G4endl
	   << "   no copies " << theNCopies << G4endl
           << "   offset1 " << theOffset1 << G4endl
           << "   offset2 " << theOffset2 << G4endl
           << "   step1 " << theStep1 << G4endl
           << "   step2 " << theStep2 << G4endl;
  }
#endif

  G4int copyNo1 = copyNo%theNCopies1;
  G4int copyNo2 = G4int(copyNo/theNCopies1);
  G4double posi1 = copyNo1*theStep1;
  G4double posi2 = copyNo2*theStep2;
  G4ThreeVector origin = posi1*theDirection1+ posi2*theDirection2;
  origin += theTranslation;

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 3 )
  {
    G4cout << " G4tgbPlaceParamSquare::ComputeTransformation() - "
           << copyNo << " = " << copyNo1 << ", X " << copyNo2 << G4endl
           << " pos: " << origin << ", axis: " << theAxis << G4endl;
  }
#endif
  //----- Set traslation and rotation
  physVol->SetTranslation(origin);
  physVol->SetCopyNo( copyNo );
  physVol->SetRotation( theRotationMatrix );
}

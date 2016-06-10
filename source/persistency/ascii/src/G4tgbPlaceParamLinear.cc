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
// $Id: G4tgbPlaceParamLinear.cc 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgbPlaceParamLinear

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgbPlaceParamLinear.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4tgrMessenger.hh"
#include "G4tgrPlaceParameterisation.hh"

// -------------------------------------------------------------------------
G4tgbPlaceParamLinear::~G4tgbPlaceParamLinear()
{
}


// -------------------------------------------------------------------------
G4tgbPlaceParamLinear::
G4tgbPlaceParamLinear( G4tgrPlaceParameterisation* tgrParam ) 
  : G4tgbPlaceParameterisation(tgrParam)
{
  //---- Get translation and rotation 
  if( tgrParam->GetParamType() == "LINEAR" )
  {
    CheckNExtraData( tgrParam, 6, WLSIZE_EQ, "G4tgbPlaceParamLinear:");
    theDirection = G4ThreeVector( tgrParam->GetExtraData()[3],
                                  tgrParam->GetExtraData()[4],
                                  tgrParam->GetExtraData()[5] );
    theAxis = kZAxis;
  }
  else
  {
    CheckNExtraData( tgrParam, 3, WLSIZE_EQ, "G4tgbPlaceParamLinear:");
    if( tgrParam->GetParamType() == "LINEAR_X" ) {
      theDirection = G4ThreeVector(1.,0.,0.);
      theAxis = kXAxis;
    } else if( tgrParam->GetParamType() == "LINEAR_Y" ) {
      theDirection = G4ThreeVector(0.,1.,0.);
      theAxis = kYAxis;
    } else if( tgrParam->GetParamType() == "LINEAR_Z" ) {
      theDirection = G4ThreeVector(0.,0.,1.);
      theAxis = kZAxis;
    }
  }

  if( theDirection.mag() == 0. )
  {
    G4Exception("G4tgbPlaceParamLinear::G4tgbPlaceParamLinear()",
                "InvalidSetup", FatalException, "Direction is zero !");
  }
  else
  {
    theDirection /= theDirection.mag();
  }

  theNCopies = G4int(tgrParam->GetExtraData()[0]);
  theStep = tgrParam->GetExtraData()[1];
  theOffset = tgrParam->GetExtraData()[2];

  theTranslation = G4ThreeVector(0.,0.,0.)+theOffset*theDirection;
    
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << " G4tgbPlaceParamLinear::G4tgbPlaceParamLinear(): "
           << " param type " << tgrParam->GetParamType() << G4endl
           << "   N copies " << theNCopies << G4endl
           << "   step " << theStep << G4endl
           << "   offset " << theOffset  << G4endl
           << "   translation " << theTranslation << G4endl
           << "   direction " << theDirection << G4endl
           << "   axis " << theAxis << G4endl;
  } 
#endif
}


// -------------------------------------------------------------------------
void G4tgbPlaceParamLinear::
ComputeTransformation(const G4int copyNo, G4VPhysicalVolume *physVol) const
{
  G4ThreeVector origin = theTranslation + copyNo*theStep*theDirection;

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 3 )
  { 
    G4cout << " G4tgbPlaceParamLinear::ComputeTransformation() -"
	   << physVol->GetName() << G4endl
           << " copyNo " << copyNo << " pos " << origin << G4endl;
  }
#endif
  //----- Set traslation and rotation
  physVol->SetTranslation(origin);
  physVol->SetCopyNo( copyNo );
  physVol->SetRotation( theRotationMatrix );
}

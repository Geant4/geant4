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
/// \file optical/wls/src/WLSUserTrackInformation.cc
/// \brief Implementation of the WLSUserTrackInformation class
//
//
//

#include "G4ios.hh"
#include "G4ThreeVector.hh"

#include "WLSUserTrackInformation.hh"

WLSUserTrackInformation::WLSUserTrackInformation ()
{
   status = undefined;
   exitPosition = G4ThreeVector(0.,0.,0.);
}

WLSUserTrackInformation::~WLSUserTrackInformation () { }

// Try adding a status flag and return if it is successful or not
// Cannot Add Undefine or a flag that conflicts with another flag
// Return true if the addition of flag is successful, false otherwise

G4bool WLSUserTrackInformation::AddStatusFlag(TrackStatus s)
{
   switch (s) {

      case left:
      case right:

        // Allow the user to set left or right
        // only if the track is undefined
        if (status == undefined) return status |= s;

        return false;
 
      case EscapedFromSide:
      case EscapedFromReadOut:

        // Allow the user to set escaped flag
        // only if the photon hasn't exited the fiber yet

        if ((status == undefined) || (status & OutsideOfFiber)) return false;

        return status |= s;
 
      case ReflectedAtMirror:
      case ReflectedAtReadOut:
      case murderee:

        return status |= s;

      case InsideOfFiber:
 
        return ( status =
         (status & ~(EscapedFromSide + EscapedFromReadOut + OutsideOfFiber)) | s);

      case OutsideOfFiber:

        return ( status = (status & ~InsideOfFiber) | s );
 
      default:
 
        return false;
   }
}

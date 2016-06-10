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
// $Id: WLSUserTrackInformation.cc 69561 2013-05-08 12:25:56Z gcosmo $
//
/// \file optical/wls/src/WLSUserTrackInformation.cc
/// \brief Implementation of the WLSUserTrackInformation class
//
//
#include "G4ios.hh"
#include "G4ThreeVector.hh"

#include "WLSUserTrackInformation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSUserTrackInformation::WLSUserTrackInformation ()
{
   fStatus = undefined;
   fExitPosition = G4ThreeVector(0.,0.,0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSUserTrackInformation::~WLSUserTrackInformation () { }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool WLSUserTrackInformation::AddStatusFlag(TrackStatus s)
// Try adding a status flag and return if it is successful or not
// Cannot Add Undefine or a flag that conflicts with another flag
// Return true if the addition of flag is successful, false otherwise
{
   switch (s) {

      case left:
      case right:

        // Allow the user to set left or right
        // only if the track is undefined
        if (fStatus == undefined) return fStatus |= s;

        return false;
 
      case EscapedFromSide:
      case EscapedFromReadOut:

        // Allow the user to set escaped flag
        // only if the photon hasn't exited the fiber yet

        if ((fStatus == undefined) || (fStatus & OutsideOfFiber)) return false;

        return fStatus |= s;
 
      case ReflectedAtMirror:
      case ReflectedAtReadOut:
      case murderee:

        return fStatus |= s;

      case InsideOfFiber:
 
        return ( fStatus =
         (fStatus & ~(EscapedFromSide + EscapedFromReadOut + OutsideOfFiber)) | s);

      case OutsideOfFiber:

        return ( fStatus = (fStatus & ~InsideOfFiber) | s );
 
      default:
 
        return false;
   }
}

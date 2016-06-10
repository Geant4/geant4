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
// $Id: WLSUserTrackInformation.hh 69561 2013-05-08 12:25:56Z gcosmo $
//
/// \file optical/wls/include/WLSUserTrackInformation.hh
/// \brief Definition of the WLSUserTrackInformation class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef WLSUserTrackInformation_h
#define WLSUserTrackInformation_h 1

#include "G4VUserTrackInformation.hh"

#include "G4ThreeVector.hh"

enum TrackStatus { undefined=0,
                   left=1, right=2, defined = 3,
                   EscapedFromSide=4, EscapedFromReadOut=8,
                   ReflectedAtMirror=16, ReflectedAtReadOut=32,
                   murderee=64, InsideOfFiber=128, OutsideOfFiber=256};

/*TrackStatus:
  undefined:
  left:                   track is going -z
  right:                  track is going +z
  defined:                left or right flag is on (Can't be Set)
  EscapedFromSide:        photon escaped through the side of the fiber
  EscapedFromReadOut:     photon escaped through the read-out end of the fiber
  ReflectedAtMirror:      photon has been reflected by the mirror at far end
  ReflectedAtReadOut:     photon has been reflected at the read-out end
  murderee                photon is artificially killed
  InsideOfFiber           Flag is on if the photon is inside of fiber
  OutsideOfFiber          Flag is on if the photon is outside of fiber
*/

class WLSUserTrackInformation : public G4VUserTrackInformation
{

  public:

    WLSUserTrackInformation();
    virtual ~WLSUserTrackInformation();
 
    const G4ThreeVector& GetExitPosition() const { return fExitPosition; }
    void SetExitPosition (const G4ThreeVector& pos) { fExitPosition = pos; }

    // Try adding a status flag and return if it is successful or not
    // Cannot Add Undefine or a flag that conflicts with another flag
    // Return true if the addition of flag is successful, false otherwise
    G4bool AddStatusFlag(TrackStatus s);

    // Check if a certain flag is on
    G4bool isStatus(TrackStatus s)
       { return s == undefined ? !(fStatus &= defined) : fStatus & s; }

  private:

    G4int fStatus;
    G4ThreeVector fExitPosition;

};

#endif

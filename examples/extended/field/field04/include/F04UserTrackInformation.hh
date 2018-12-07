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
/// \file field/field04/include/F04UserTrackInformation.hh
/// \brief Definition of the F04UserTrackInformation class
//

#ifndef F04UserTrackInformation_h
#define F04UserTrackInformation_h 1

#include "G4VUserTrackInformation.hh"

enum TrackStatus { undefined, left, right, reverse };

/*TrackStatus:
  undefined:
  left:      track is going -z
  right:     track is going +z
  reverse:   track has reversed direction from left to right
*/

class F04UserTrackInformation : public G4VUserTrackInformation
{

  public:

    F04UserTrackInformation();
    virtual ~F04UserTrackInformation();

    void SetTrackStatusFlag(TrackStatus s){ fStatus = s; }
    TrackStatus GetTrackStatusFlag()const { return fStatus; }

  private:

    TrackStatus fStatus;

};

#endif

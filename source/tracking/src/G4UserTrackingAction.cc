// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UserTrackingAction.cc,v 1.5 2001-02-08 07:39:53 tsasaki Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// ---------------------------------------------------------------
//
// G4UserTrackingAction.cc
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
// ---------------------------------------------------------------

#include "G4UserTrackingAction.hh"

/////////////////////////////////////////////////////////
void G4UserTrackingAction::
     SetTrackingManagerPointer(G4TrackingManager* pValue)
/////////////////////////////////////////////////////////
{
  fpTrackingManager = pValue;
}  









// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TrackStatus.hh,v 1.1 1999-01-07 16:14:24 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// G4TrackStatus.hh
//
// Description:
//   This is an enumerator to define the current status
//   of the track which is under the transportation.
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

#ifndef G4TrackStatus_h
#define G4TrackStatus_h 1

//////////////////
enum G4TrackStatus
//////////////////
{

  fAlive,             // Continue the tracking
  fStopButAlive,      // Invoke active rest physics processes and
                      // and kill the current track afterward
  fStopAndKill,       // Kill the current track

  fKillTrackAndSecondaries,
                      // Kill the current track and also associated
                      // secondaries.
  fSuspend,           // Suspend the current track
  fPostponeToNextEvent
                      // Postpones the tracking of thecurrent track 
                      // to the next event.

};

#endif



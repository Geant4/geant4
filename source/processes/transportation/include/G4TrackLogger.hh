//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4TrackLogger.hh,v 1.1 2002-07-10 15:51:04 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4TrackLogger
//
// Class description:
// 
// This class loggs every track via it's id. It may tell if
// a track has been logged. The loggs are cleared when ever 
// the object has been informed about a new event.
// 
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
//

#ifndef G4TrackLogger_hh
#define G4TrackLogger_hh G4TrackLogger_hh

#include "globals.hh"
#include "g4std/set"

class G4TrackLogger {
public:
  G4TrackLogger();
  ~G4TrackLogger();

  void SetEventID(G4int id);
    // inform the object about the event number
    // if the event number changes the loggs are cleared.

  G4bool FirstEnterance(G4int trid);
    // returns true if the track is new to this event.

  typedef G4std::set<G4int > TrackIDsSet;
    // log for the track ids.

private:
  G4int fPreviousEventID;
  TrackIDsSet fTrackIDsSet;

};

#endif

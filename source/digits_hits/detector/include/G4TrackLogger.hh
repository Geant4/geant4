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
// $Id: G4TrackLogger.hh 67992 2013-03-13 10:59:57Z gcosmo $
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
#include <set>

class G4TrackLogger {
public:
  G4TrackLogger();
  ~G4TrackLogger();

  void SetEventID(G4int id);
    // inform the object about the event number
    // if the event number changes the loggs are cleared.

  G4bool FirstEnterance(G4int trid);
    // returns true if the track is new to this event.

  typedef std::set<G4int > TrackIDsSet;
    // log for the track ids.

private:
  G4int fPreviousEventID;
  TrackIDsSet fTrackIDsSet;

};

#endif

// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20TrackingAction.hh,v 1.1 2001-05-25 12:50:06 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//-----------------------------------------------------------------------------

#ifndef Tst20TrackingAction_h
#define Tst20TrackingAction_h 1

#include "G4UserTrackingAction.hh"

class Tst20TrackingAction : public G4UserTrackingAction {

  public:
    Tst20TrackingAction(){};
    virtual ~Tst20TrackingAction(){};
   
    virtual void PreUserTrackingAction(const G4Track* aTrack);

};

#endif

// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst14TrackingAction.hh,v 1.2 2001-05-23 11:01:21 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//-----------------------------------------------------------------------------

#ifndef Tst14TrackingAction_h
#define Tst14TrackingAction_h 1

#include "G4UserTrackingAction.hh"

class Tst14TrackingAction : public G4UserTrackingAction {

  public:
    Tst14TrackingAction(){};
    virtual ~Tst14TrackingAction(){};
   
    virtual void PreUserTrackingAction(const G4Track* aTrack);

};

#endif

// $Id: A01TrackingAction.hh,v 1.1 2002-11-13 07:19:00 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#ifndef A01TrackingAction_h
#define A01TrackingAction_h 1

#include "G4UserTrackingAction.hh"

class A01TrackingAction : public G4UserTrackingAction
{
  public:
    A01TrackingAction();
    virtual ~A01TrackingAction();

    virtual void PreUserTrackingAction(const G4Track*);
};

#endif

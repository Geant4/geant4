// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ClassificationOfNewTrack.hh,v 1.3 1999-12-15 14:49:38 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#ifndef G4ClassificationOfNewTrack_hh
#define G4ClassificationOfNewTrack_hh 1

// class description:
//
//  This header file contain an enumeration for the possible classifications
// for trackes newly pushed to the stack. G4UserStackingAction can set the
// classification.

enum G4ClassificationOfNewTrack
{ 
  fUrgent,    // put into the urgent stack
  fWaiting,   // put into the waiting stack
  fPostpone,  // postpone to the next event
  fKill       // kill
};

#endif


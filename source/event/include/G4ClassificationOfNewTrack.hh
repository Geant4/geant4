// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ClassificationOfNewTrack.hh,v 1.1 1999-01-07 16:06:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#ifndef G4ClassificationOfNewTrack_hh
#define G4ClassificationOfNewTrack_hh 1

enum G4ClassificationOfNewTrack
{ 
  fUrgent,    // put into the urgent stack
  fWaiting,   // put into the waiting stack
  fPostpone,  // postpone to the next event
  fKill       // kill
};

#endif


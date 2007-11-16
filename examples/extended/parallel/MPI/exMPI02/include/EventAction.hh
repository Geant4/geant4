// $Id: EventAction.hh,v 1.1 2007-11-16 14:29:32 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   EventAction.hh
//
//                                         2007 Q
// ====================================================================
#ifndef EVENT_ACTION_H
#define EVENT_ACTION_H

#include "G4UserEventAction.hh"

// ====================================================================
//
// class definition
//
// ====================================================================

class EventAction : public G4UserEventAction {

public:
  EventAction();
  ~EventAction();

  virtual void BeginOfEventAction(const G4Event* aevent);
  virtual void EndOfEventAction(const G4Event* aevent);

};

#endif

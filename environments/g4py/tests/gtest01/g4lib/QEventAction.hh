// $Id: QEventAction.hh,v 1.1 2006-05-11 03:00:08 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   QEventAction.hh
//
//                                         2005 Q
// ====================================================================
#ifndef Q_EVENT_ACTION_H
#define Q_EVENT_ACTION_H

#include "G4UserEventAction.hh"

// ====================================================================
//
// class definition
//
// ====================================================================

class QEventAction : public G4UserEventAction {
public:
  QEventAction();
  ~QEventAction();

  virtual void BeginOfEventAction(const G4Event* aevent);
  virtual void EndOfEventAction(const G4Event* aevent);

};

#endif

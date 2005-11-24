//
// File name:     RadmonApplicationEventNumbering.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationEventNumbering.cc,v 1.1 2005-11-24 02:34:21 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonApplicationEventNumbering.hh"
#include "G4Event.hh"

void                                            RadmonApplicationEventNumbering :: OnBeginOfEvent(const G4Event * event)
{
 if (dumpEvery==0)
  return;
  
 G4int eventId(event->GetEventID());

 if ((eventId%dumpEvery) == 0)
  G4cout << "RadmonApplicationEventNumbering::OnBeginOfEvent: Begin of event " << eventId << '.' << G4endl;
}

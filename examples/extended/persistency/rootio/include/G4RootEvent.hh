// $Id: G4RootEvent.hh,v 1.1 2002-12-04 02:44:28 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4RootEvent.hh
//
// History:
//   '02.5.7  Youhei Morita  Initial creation

#ifndef ROOT_FADS_EVENT_HH
#define ROOT_FADS_EVENT_HH 1

class G4Event;
class GenEvent;
class G4MCTEvent;
class G4Pevent;

// Class inherited:
#include "TObject.h"

// Class Description:
//   FADS persistent event object.
// 
//   This class has pointers to HepMCGenEvent, MCTruth and G4Event.

class G4RootEvent
 : public TObject
{
    public: // With description
      G4RootEvent() {};
      // Constructor

      virtual ~G4RootEvent();
      // Destructor

    public: // With description
      int GetEventID() { return m_id; };
      // returns the event ID.

      G4Pevent* MakeTransientObject();
      // Make transient G4Pevent from G4RootEvent

    private:
      int m_id;

    ClassDef (G4RootEvent,1)

}; // End of class G4RootEvent

#endif


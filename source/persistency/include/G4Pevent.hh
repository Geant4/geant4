// $Id: G4Pevent.hh,v 1.1 2002-11-24 13:45:23 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4Pevent.hh
//
// History:
//   '01.11.18  Youhei Morita  Initial creation

#ifndef G4PEVENT_HH
#define G4PEVENT_HH 1

#include "G4Event.hh"
// not yet // #include "G4MCTEvent.hh"
#include "CLHEP/HepMC/GenEvent.h"

// Class Description:
//   Geant4 event object for store.
// 
//   This class has pointers to HepMC::GenEvent, MCTruth and G4Event.
// 
//   In the event store operation, this object will be created in the concrete
//   class of G4VPEventIO::Store() method, and will be deleted immediately after
//   creating persistent Geant4 event. 
// 
//   In the event retrieve operation, this object will be created in the
//   concrete Persistency::Retrieve() method.  The retrieved G4Pevent
//   has to be deleted by the user method which called the Retrieve().

class G4Pevent
{
    public: // With description
      // G4Pevent( HepMC::GenEvent* hepevt, G4MCTEvent* mctevt, G4Event* g4evt );
      G4Pevent( HepMC::GenEvent* hepevt, G4Event* g4evt );
      // Constructor

      ~G4Pevent();
      // Destructor

    public: // With description
      int GetEventID() { return m_id; };
      // returns the event ID.

      G4Event* GetG4Event() { return f_g4evt; };
      // returns the G4Event.

      HepMC::GenEvent* GetHepMCGenEvent() { return f_hepevt; };
      // returns the HepMC GenEvent.

      // G4MCTEvent* GetG4MCTEvent() { return f_mctevt; };
      // returns the MCTruth event.

      int GetGenEventID() const { return genEventID; };
      // returns the GenEvent ID.

      void SetGenEventID(int id) { genEventID=id; };
      // set the GenEvent ID.

    private:
      HepMC::GenEvent* f_hepevt;
      int       genEventID;
      // G4MCTEvent* f_mctevt;
      G4Event*  f_g4evt;
      int       m_id;

}; // End of class G4Pevent

#endif


// $Id: G4VPEventIO.hh,v 1.4 2002-12-04 11:44:06 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4VPEventIO.hh
//
// History:
//   '01.08.10  Youhei Morita  Initial creation (with "fadsclass3")

#ifndef V_P_EVENT_I_O_HH
#define V_P_EVENT_I_O_HH 1

#include "G4Event.hh"

#ifndef WIN32
  #include "CLHEP/HepMC/GenEvent.h"
  // not yet // #include "G4MCTEvent.hh"
#include "G4Pevent.hh"
#endif

// Class Description:
//   Abstract base class for storing and retrieving an event object

class G4VPEventIO
{
    public: // With description
      G4VPEventIO();
      // Constructor

      virtual ~G4VPEventIO() {};
      // Destructor

    public: // With description
      void SetVerboseLevel(int v) { m_verbose = v; };
      // Set verbose level.

      inline G4int CurrentEventID() { return m_currentEvtID; };
      // Returns the current event id.

      virtual G4bool Store( const G4Event* anEvent ) =0;
      // Store a Geant4 event.
#ifndef WIN32
      // virtual G4bool Store( HepMC::GenEvent* hepevt, G4MCTEvent* mctevt, const G4Event* anEvent ) =0;
      virtual G4bool Store( HepMC::GenEvent* hepevt, const G4Event* anEvent ) =0;
      // Store a Geant4 event.

      virtual G4bool Retrieve( G4Pevent*& anEvent ) =0;
      // Retrieve a Geant4 event.
#endif
      virtual G4bool Retrieve( G4Event*& anEvent ) =0;
      // Retrieve a Geant4 event.

    protected:
      G4int m_verbose;
      G4int m_currentEvtID;

}; // End of class G4VPEventIO

#endif


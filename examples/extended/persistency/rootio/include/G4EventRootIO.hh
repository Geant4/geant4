// $Id: G4EventRootIO.hh,v 1.1 2002-12-04 02:44:28 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4EventRootIO.hh
//
// History:
//   '02.5.7  Youhei Morita  Initial creation

#ifndef EVENT_ROOT_IO_HH
#define EVENT_ROOT_IO_HH 1

#include "G4Event.hh"
#include "G4PersistencyCenter.hh"
#include "G4RootIOManager.hh"
#include "G4HitRootIO.hh"
#include "G4DigitRootIO.hh"
#include "G4HepMCRootIO.hh"
#include "G4MCTruthRootIO.hh"
#include "G4RootEvent.hh"

// Class inherited:
#include "G4VPEventIO.hh"

// Class Description:
//   Manager class to store and retrieve G4Event objects

class G4EventRootIO
 : public G4VPEventIO
{
    public: // With description
      G4EventRootIO( G4PersistencyManager* pc );
      // Constructor

      virtual ~G4EventRootIO(){};
      // Destructor

    public: // With description
      // G4bool Store( HepMC::GenEvent* hepevt, G4MCTEvent* mctevt, const G4Event* g4evt);
      G4bool Store( HepMC::GenEvent* hepevt, const G4Event* g4evt);
      // Store a FADS event.

      G4bool Store( const G4Event* anEvent );
      // Store a Geant4 event.

      G4bool Retrieve( G4Pevent*& anEvent );
      // Retrieve a FADS event.

      G4bool Retrieve( G4Event*& anEvent );
      // Retrieve a Geant4 event.

    private:
      G4PersistencyManager* f_pc;
      G4HepMCRootIO*        f_hepmcio;
      G4MCTruthRootIO*      f_mctruthio;
      G4HitRootIO*          f_PHCMan;
      G4DigitRootIO*        f_PDCMan;

}; // End of class G4EventRootIO

#endif


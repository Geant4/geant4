// $Id: G4HepMCRootIO.hh,v 1.1 2002-12-04 02:44:28 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4HepMCRootIO.hh
//
// History:
//   '02.5.7  Youhei Morita  Initial creation

#ifndef HEPMC_ROOT_IO_HH
#define HEPMC_ROOT_IO_HH 1

#include "CLHEP/HepMC/GenEvent.h"
#include "CLHEP/HepMC/GenParticle.h"
#include "G4PersistencyCenter.hh"
#include "G4RootIOManager.hh"

// Class inherited:
#include "G4VHepMCIO.hh"

// Class Description:
//   Concreate class for storing and retrieving HepMC events.

class G4HepMCRootIO
 : public G4VHepMCIO
{
    public: // With description
      G4HepMCRootIO() {};
      // Constructor

      virtual ~G4HepMCRootIO() {};
      // Destructor

    public: // With description
      bool Store(HepMC::GenEvent* evt);
      // Method for storing HepMC GenEvent.

      bool Retrieve(HepMC::GenEvent*& evt, int id=-1);
      // Method for retrieving HepMC GenEvent.

      int LastEventID() { return m_lastid; };
      // returns the last event ID.

    private:
      HepMC::GenEvent*       f_lastevt;
      int                    m_lastid;

}; // End of class G4HepMCRootIO

#endif


// $Id: G4VHepMCIO.hh,v 1.2 2002-12-04 10:25:49 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4VHepMCIO.hh
//
// History:
//   '01.11.18  Youhei Morita  Initial creation

#ifndef G4V_HEPMC_IO_HH
#define G4V_HEPMC_IO_HH 1

#include "G4Types.hh"
#include "CLHEP/HepMC/GenEvent.h"

// Class Description:
//   Abstract base class for storing and retrieving HepMC events.

class G4VHepMCIO
{
    public: // With description
      G4VHepMCIO();
      // Constructor

      virtual ~G4VHepMCIO() {};
      // Destructor

    public: // With description
      virtual G4bool Store(HepMC::GenEvent*) =0;
      // Pure virtual method for storing HepMC GenEvent.
      // Each persistency package should implement a concrete method
      // of storing the HepMC::GenEvent with this signature.

      virtual G4bool Retrieve(HepMC::GenEvent*&, int id=-1) =0;
      // Pure virtual method for retrieving HepMC GenEvent.
      // Each persistency package should implement a concrete method
      // of storing the HepMC::GenEvent with this signature.

      virtual int LastEventID() =0;
      //  Pure virtual method for the event ID

      void SetVerboseLevel(int v) { m_verbose = v; };
      // Set verbose level.

    protected:
      G4int m_verbose;

}; // End of class G4VHepMCIO

#endif


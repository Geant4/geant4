//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
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


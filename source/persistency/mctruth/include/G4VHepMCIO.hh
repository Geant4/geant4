//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
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


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
// INCL++ intra-nuclear cascade model
// Alain Boudard, CEA-Saclay, France
// Joseph Cugnon, University of Liege, Belgium
// Jean-Christophe David, CEA-Saclay, France
// Pekka Kaitaniemi, CEA-Saclay, France, and Helsinki Institute of Physics, Finland
// Sylvie Leray, CEA-Saclay, France
// Davide Mancusi, CEA-Saclay, France
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#ifndef G4INCLPauliStrictStandard_hh
#define G4INCLPauliStrictStandard_hh 1

#include "G4INCLIPauli.hh"
#include "G4INCLNucleus.hh"

namespace G4INCL {

  class PauliStrictStandard : public IPauli {
  public:
    PauliStrictStandard();
    ~PauliStrictStandard();

    /// \brief Dummy copy constructor to silence Coverity warning
    PauliStrictStandard(const PauliStrictStandard &rhs);

    /// \brief Dummy assignment operator to silence Coverity warning
    PauliStrictStandard &operator=(const PauliStrictStandard &rhs);

    G4bool isBlocked(ParticleList const &, Nucleus const * const);

  private:
    IPauli *theStrictBlocker, *theStandardBlocker;
  };
}

#endif

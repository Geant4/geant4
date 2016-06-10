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

#include "G4INCLCDPP.hh"
#include <functional>
#include <algorithm>

namespace G4INCL {

  CDPP::CDPP() :
    Sk(0.0),
    TbelowTf(0.0),
    thePotential(NULL)
  {}

  CDPP::~CDPP() {}

  G4bool CDPP::isBlocked(ParticleList const &created, Nucleus const * const nucleus) {
    G4double S = nucleus->computeSeparationEnergyBalance();

    thePotential = nucleus->getPotential();

    ParticleList const &remnantParticles = nucleus->getStore()->getParticles();

    Sk = 0.0;
    TbelowTf = 0.0;

    std::for_each(remnantParticles.begin(), remnantParticles.end(), std::bind1st(std::mem_fun(&G4INCL::CDPP::processOneParticle), this));
    std::for_each(created.begin(), created.end(), std::bind1st(std::mem_fun(&G4INCL::CDPP::processOneParticle), this));

    const G4double Tinitial = nucleus->getInitialInternalEnergy();
    const G4double Eblock = TbelowTf - Tinitial - Sk - S;

    return (Eblock < 0.0);
  }
}

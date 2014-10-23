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

#include "G4INCLPauliGlobal.hh"
#include "G4INCLRandom.hh"

namespace G4INCL {

  PauliGlobal::PauliGlobal() {}
  PauliGlobal::~PauliGlobal() {}

  G4bool PauliGlobal::isBlocked(ParticleList const &pL, Nucleus const * const n) {
    for(ParticleIter p=pL.begin(), e=pL.end(); p!=e; ++p) {
      // Pauli blocking only applies to nucleons
      if(!(*p)->isNucleon()) continue;

      // If the particle is above T_F, it is never blocked
      const ParticleType t = (*p)->getType();
      const G4double pFermi = n->getPotential()->getFermiMomentum(t);
      const G4double pFermiSquared = pFermi*pFermi;
      if((*p)->getMomentum().mag2() > pFermiSquared) continue;

      // Count particles of the same type as p below the Fermi sea
      ParticleList const &particles = n->getStore()->getParticles();
      G4int nSea = 0;
      for(ParticleIter i=particles.begin(), end=particles.end(); i!=end; ++i) {
        if((*i)->getType() != t) continue;
        const G4double pmod2 = (*i)->getMomentum().mag2();
        if(pmod2<pFermiSquared) nSea++;
      }

      // Compute the blocking probability
      G4double probBlocking;
      if(t==Proton)
        probBlocking = ((G4double) nSea)/((G4double) n->getInitialZ());
      else
        probBlocking = ((G4double) nSea)/((G4double) (n->getInitialA() - n->getInitialZ()));

      // The avatar is blocked if any particle is blocked
      if(Random::shoot() < probBlocking) return true;

    }

    // Not blocked
    return false;

  }
}

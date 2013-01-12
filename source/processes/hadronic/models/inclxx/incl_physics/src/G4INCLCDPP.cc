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
// Pekka Kaitaniemi, CEA and Helsinki Institute of Physics
// Davide Mancusi, CEA
// Alain Boudard, CEA
// Sylvie Leray, CEA
// Joseph Cugnon, University of Liege
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4INCLCDPP.hh"

namespace G4INCL {

  CDPP::CDPP() {}

  CDPP::~CDPP() {}

  G4bool CDPP::isBlocked(ParticleList const created, Nucleus const * const nucleus) const {
    G4double S = nucleus->computeSeparationEnergyBalance();

    const G4double Sp = nucleus->getPotential()->getSeparationEnergy(Proton);
    const G4double Sn = nucleus->getPotential()->getSeparationEnergy(Neutron);

    ParticleList remnantParticles = nucleus->getStore()->getParticles();
    remnantParticles.insert(remnantParticles.end(), created.begin(), created.end());

    G4double Sk = 0.0;
    G4double TbelowTf = 0.0;
    for(ParticleIter i = remnantParticles.begin(); i != remnantParticles.end();
        ++i) {

      if((*i)->isNucleon()) {
        const G4double Tf = nucleus->getPotential()->getFermiEnergy(*i);
        const G4double T = (*i)->getKineticEnergy();

        if(T > Tf) {
          const G4double sep = nucleus->getPotential()->getSeparationEnergy(*i);
          Sk += sep;
        } else {
          TbelowTf += T - (*i)->getPotentialEnergy();
        }
      } else if((*i)->isResonance()) {
        const G4double Tf = nucleus->getPotential()->getFermiEnergy(*i);
        const G4double T = (*i)->getKineticEnergy();

        if(T > Tf) {
          const G4double sep = nucleus->getPotential()->getSeparationEnergy(*i);
          Sk += sep;
        } else { // Ugly! We should use total energies everywhere!
          TbelowTf += (*i)->getEnergy() - ParticleTable::getINCLMass(Proton) - (*i)->getPotentialEnergy();
        }
      } else if((*i)->getType() == PiPlus)
        Sk += Sp - Sn;
      else if((*i)->getType() == PiMinus)
        Sk += Sn - Sp;

    }
    G4double Tinitial = nucleus->getInitialInternalEnergy();
    G4double Eblock = TbelowTf - Tinitial - Sk - S;

    return (Eblock < 0.0);
  }
}

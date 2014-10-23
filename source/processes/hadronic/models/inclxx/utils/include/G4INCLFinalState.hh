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

/*
 * G4INCLIChannel.h
 *
 *  \date Jun 5, 2009
 * \author Pekka Kaitaniemi
 */

#ifndef G4INCLFINALSTATE_H_
#define G4INCLFINALSTATE_H_

#include "G4INCLParticle.hh"
#include <string>

namespace G4INCL {

  enum FinalStateValidity {
    ValidFS,
    PauliBlockedFS,
    NoEnergyConservationFS,
    ParticleBelowFermiFS,
    ParticleBelowZeroFS
  };

  /**
   * Final state of an interaction
   */
  class FinalState {
  public:
    FinalState();
    virtual ~FinalState();

    void reset();

    void setTotalEnergyBeforeInteraction(G4double E) { totalEnergyBeforeInteraction = E; };
    G4double getTotalEnergyBeforeInteraction() const { return totalEnergyBeforeInteraction; };

    void addModifiedParticle(Particle *p);
    void addOutgoingParticle(Particle *p);
    void addDestroyedParticle(Particle *p);
    void addCreatedParticle(Particle *p);
    void addEnteringParticle(Particle *p);

    ParticleList const &getModifiedParticles() const;
    ParticleList const &getOutgoingParticles() const;
    ParticleList const &getDestroyedParticles() const;
    ParticleList const &getCreatedParticles() const;
    ParticleList const &getEnteringParticles() const;

    FinalStateValidity getValidity() const { return validity; }
    void makeValid() { validity = ValidFS; }
    void makePauliBlocked() { validity = PauliBlockedFS; }
    void makeNoEnergyConservation() { validity = NoEnergyConservationFS; }
    void makeParticleBelowFermi() { validity = ParticleBelowFermiFS; }
    void makeParticleBelowZero() { validity = ParticleBelowZeroFS; }

    std::string print() const;

  private:
    ParticleList outgoing, created, destroyed, modified, entering;
    G4double totalEnergyBeforeInteraction;
    FinalStateValidity validity;
  };

}

#endif /* G4INCLFINALSTATE_H_ */

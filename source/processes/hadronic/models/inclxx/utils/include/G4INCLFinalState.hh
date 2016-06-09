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
// INCL++ revision: v5.0_rc3
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/*
 * G4INCLIChannel.h
 *
 *  Created on: Jun 5, 2009
 *      Author: Pekka Kaitaniemi
 */

#ifndef G4INCLFINALSTATE_H_
#define G4INCLFINALSTATE_H_

#include "G4INCLParticle.hh"
#include <string>

namespace G4INCL {

  enum FinalStateValidity {ValidFS, PauliBlockedFS, NoEnergyConservationFS};

  /**
   * Final state of an G4interaction
   */
  class FinalState {
  public:
    FinalState();
    virtual ~FinalState();

    void setTotalEnergyBeforeInteraction(G4double E) { totalEnergyBeforeInteraction = E; };
    G4double getTotalEnergyBeforeInteraction() const { return totalEnergyBeforeInteraction; };

    void addModifiedParticle(Particle *p);
    void addOutgoingParticle(Particle *p);
    void addDestroyedParticle(Particle *p);
    void addCreatedParticle(Particle *p);

    ParticleList getModifiedParticles() const;
    ParticleList getOutgoingParticles() const;
    ParticleList getDestroyedParticles() const;
    ParticleList getCreatedParticles() const;

    void setBlockedDelta(Particle * const p) { blockedDelta = p; }
    Particle *getBlockedDelta() { return blockedDelta; }

    FinalStateValidity getValidity() const { return validity; }
    void makeValid() { validity = ValidFS; }
    void makePauliBlocked() { validity = PauliBlockedFS; }
    void makeNoEnergyConservation() { validity = NoEnergyConservationFS; }

    std::string prG4int() const;

  private:
    ParticleList outgoing, created, destroyed, modified;
    G4double totalEnergyBeforeInteraction;
    FinalStateValidity validity;
    Particle *blockedDelta;
  };

}

#endif /* G4INCLFINALSTATE_H_ */

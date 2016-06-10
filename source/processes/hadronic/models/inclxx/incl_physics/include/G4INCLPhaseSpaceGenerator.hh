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

#ifndef G4INCLPHASESPACEGENERATOR_HH
#define G4INCLPHASESPACEGENERATOR_HH 1

#include "G4INCLIPhaseSpaceGenerator.hh"
#include "G4INCLConfig.hh"

namespace G4INCL {
  namespace PhaseSpaceGenerator {
    /// \brief Generate an event in the CM system
    void generate(const G4double sqrtS, ParticleList &particles);

    /** \brief Generate a biased event in the CM system
     *
     * This method first generates a "flat" event by calling generate(). The
     * particles are subsequently rotated in such a way that one of them
     * (identified by the parameter index) is biased towards the collisionAxis
     * with an exponential distribution of the form
     * \f[
     * \exp(B\cdot t)
     * \f]
     * where \f$t\f$ is the usual Mandelstam variable. The incoming momentum
     * is taken to be the momentum of particles[index] at the moment of the
     * call.
     *
     * \param sqrtS total energy in the centre of mass, in MeV
     * \param particles list of particles for which the event will be generated
     *        (modified on exit)
     * \param index index of the particle to be biased; all the other particles
     *        will follow
     * \param slope slope \f$B\f$ of the angular distribution: \f$\exp(Bt)\f$,
     *        in (GeV/c)^(-2)
     */
    void generateBiased(const G4double sqrtS, ParticleList &particles, const size_t index, const G4double slope);

    void setPhaseSpaceGenerator(IPhaseSpaceGenerator *g);

    IPhaseSpaceGenerator *getPhaseSpaceGenerator();

    void deletePhaseSpaceGenerator();

    void initialize(Config const * const theConfig);
  }
}

#endif

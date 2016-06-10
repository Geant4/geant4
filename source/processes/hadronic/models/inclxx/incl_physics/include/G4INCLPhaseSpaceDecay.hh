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

#ifndef G4INCLPHASESPACEDECAY_HH_
#define G4INCLPHASESPACEDECAY_HH_

#include "G4INCLThreeVector.hh"
#include "G4INCLParticle.hh"
#include <vector>

namespace G4INCL {
  namespace PhaseSpaceDecay {

    /** \brief Generate decay momenta according to a uniform phase-space model
     *
     * This function will assign momenta to the particles in the list that is
     * passed as an argument.
     *
     * \param initialMass mass of the decaying system
     * \param theBoostVector boost vector of the decaying system
     * \param particles list of decay particles
     */
    void decay(G4double initialMass, const ThreeVector &theBoostVector, ParticleList &particles);

  }
}

#endif // G4INCLPHASESPACEDECAY_HH_

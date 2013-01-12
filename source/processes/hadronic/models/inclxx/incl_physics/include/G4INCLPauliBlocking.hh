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

#ifndef G4INCLPauliBlocking__hh
#define G4INCLPauliBlocking__hh 1

#include "G4INCLIPauli.hh"
#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"

namespace G4INCL {

  /**
   * Pauli blocking
   *
   */
  class Pauli {
  public:
    /** \brief Check Pauli blocking.
     *
     * Note: This is a "pure" function: it doesn't retain or modify
     * any state at all and thus only depends on its arguments.
     *
     * \param p list of modified and created particles
     * \param n the nucleus
     */
    static G4bool isBlocked(ParticleList const p, Nucleus const * const n);

    /** \brief Check CDPP blocking.
     *
     * Note: This is a "pure" function: it doesn't retain or modify
     * any state at all and thus only depends on its arguments.
     *
     * \param p list of created particles
     * \param n the nucleus
     */
    static G4bool isCDPPBlocked(ParticleList const p, Nucleus const * const n);

    /**
     * Get the Pauli blocker algorithm.
     */
    static IPauli const * getBlocker() { return thePauliBlocker; }

    /**
     * Get the CDPP blocker algorithm.
     */
    static IPauli const * getCDPP() { return theCDPP; }

    /**
     * Set the Pauli blocker algorithm.
     */
    static void setBlocker(IPauli const * const);

    /**
     * Set the CDPP blocker algorithm.
     */
    static void setCDPP(IPauli const * const);

    /**
     * Delete blockers
     */
    static void deleteBlockers() {
      delete thePauliBlocker;
      thePauliBlocker=NULL;
      delete theCDPP;
      theCDPP=NULL;
    }

  protected:
    Pauli() {}
    ~Pauli() {}

  private:
    static IPauli const * thePauliBlocker;
    static IPauli const * theCDPP;
  };

}

#endif

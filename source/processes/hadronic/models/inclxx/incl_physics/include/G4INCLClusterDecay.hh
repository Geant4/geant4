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

/** \file G4INCLClusterDecay.hh
 * \brief Static class for carrying out cluster decays
 *
 * Created on: 6th July 2011
 *     Author: Davide Mancusi
 */

#ifndef G4INCLCLUSTERDECAY_HH
#define G4INCLCLUSTERDECAY_HH

#include "G4INCLCluster.hh"
#include "G4INCLParticleTable.hh"

namespace G4INCL {

  /**
   * Pauli blocking
   *
   */
  class ClusterDecay {
  public:
    /// \brief True if the cluster is stable.
    static G4bool isStable(Cluster const * const c) {
      const G4int Z = c->getZ();
      const G4int A = c->getA();
      return (ParticleTable::clusterDecayMode[Z][A]==ParticleTable::StableCluster);
    }

    /** \brief Carries out a cluster decay
     *
     * \param c cluster that should decay
     * \return list of decay products
     */
    static ParticleList decay(Cluster * const c);

  protected:
    ClusterDecay() {}
    ~ClusterDecay() {}

  private:
    /** \brief Recursively decay clusters
     *
     * \param c cluster that should decay
     * \param decayProducts decay products are appended to the end of this list
     */
    static void recursiveDecay(Cluster * const c, ParticleList *decayProducts);

    /// \brief Carries out two-body decays
    static void twoBodyDecay(Cluster * const c, ParticleTable::ClusterDecayType theDecayMode, ParticleList *decayProducts);

    /// \brief Carries out three-body decays
    static void threeBodyDecay(Cluster * const c, ParticleTable::ClusterDecayType theDecayMode, ParticleList *decayProducts);

  };

}

#endif // G4INCLCLUSTERDECAY

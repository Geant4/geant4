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

/** \file G4INCLClusterDecay.hh
 * \brief Static class for carrying out cluster decays
 *
 * \date 6th July 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLCLUSTERDECAY_HH
#define G4INCLCLUSTERDECAY_HH

#include "G4INCLCluster.hh"
#include "G4INCLParticleTable.hh"

namespace G4INCL {

  /// \brief Namespace for functions that handle decay of unstable clusters
  namespace ClusterDecay {

      /// \brief True if the cluster is stable.
      G4bool isStable(Cluster const * const c);

      /** \brief Carries out a cluster decay
       *
       * \param c cluster that should decay
       * \return list of decay products
       */
      ParticleList decay(Cluster * const c);

      enum ClusterDecayType {
        StableCluster,
        NeutronDecay,
        ProtonDecay,
        AlphaDecay,
        TwoProtonDecay,
        TwoNeutronDecay,
        ProtonUnbound,
        NeutronUnbound,
        LambdaUnbound,
        LambdaDecay
      };

      /** \brief Table for cluster decays
       *
       * Definition of "Stable": halflife > 1 ms
       *
       * These table includes decay data for clusters that INCL presently does
       * not produce. It can't hurt.
       *
       * Unphysical nuclides (A<Z) are marked as stable, but should never be
       * produced by INCL. If you find them in the output, something is fishy.
       */
      extern G4ThreadLocal ClusterDecayType clusterDecayMode[ParticleTable::clusterTableSSize][ParticleTable::clusterTableZSize][ParticleTable::clusterTableASize];

  }

}

#endif // G4INCLCLUSTERDECAY

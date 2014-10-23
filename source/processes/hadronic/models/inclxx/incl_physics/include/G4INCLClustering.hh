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

/** \file G4INCLClustering.hh
 * \brief Static class for cluster formation
 *
 * \date 13th July 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLCLUSTERING_HH
#define G4INCLCLUSTERING_HH

#include "G4INCLIClusteringModel.hh"
#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"

namespace G4INCL {

  /// \brief Cluster formation
  namespace Clustering {
    /** \brief Call the clustering algorithm
     *
     * Choose a cluster candidate to be produced. At this point we
     * don't yet decide if it can pass through the Coulomb barrier or
     * not.
     */
    Cluster* getCluster(Nucleus *n, Particle *p);

    /// \brief Determine whether the cluster can escape or not
    G4bool clusterCanEscape(Nucleus const * const n, Cluster const * const c);

    /// \brief Get the clustering model
    IClusteringModel *getClusteringModel();

    /// \brief Set the clustering model
    void setClusteringModel(IClusteringModel * const model);

    /// \brief Delete the clustering model
    void deleteClusteringModel();

    /// \brief Initialize the clustering model based on the Config object
    void initialize(Config const * const theConfig);

  }
}

#endif

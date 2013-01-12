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

  /**
   * Cluster formation
   *
   */
  class Clustering {
  public:
    /**
     * Choose a cluster candidate to be produced. At this point we
     * don't yet decide if it can pass through the Coulomb barrier or
     * not.
     */
    static Cluster* getCluster(Nucleus *n, Particle *p) {
#if !defined(NDEBUG) && !defined(INCLXX_IN_GEANT4_MODE)
      Cluster * const c=theClusteringModel->getCluster(n,p);
// assert(!c || c->getA()<=n->getA()/2);
      return c;
#else
      return theClusteringModel->getCluster(n,p);
#endif
    }

    /**
     * Determine whether cluster can escape or not.
     */
    static G4bool clusterCanEscape(Nucleus const * const n, Cluster const * const c) {
      return theClusteringModel->clusterCanEscape(n, c);
    }

    /// \brief Get the clustering model.
    static IClusteringModel *getClusteringModel() { return theClusteringModel; }

    /// \brief Set the clustering model
    static void setClusteringModel(IClusteringModel * const model) {
      theClusteringModel = model;
    }

    /**
     * Delete clustering model
     */
    static void deleteClusteringModel() {
      delete theClusteringModel;
      theClusteringModel = 0;
    }

  protected:
    Clustering() {}
    ~Clustering() {}

  private:
    static IClusteringModel *theClusteringModel;
  };

}

#endif

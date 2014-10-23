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

/** \file G4INCLClustering.cc
 * \brief Static class for cluster formation
 *
 * \date 13th July 2011
 * \author Davide Mancusi
 */

#include "G4INCLIClusteringModel.hh"
#include "G4INCLClustering.hh"
#include "G4INCLClusteringModelIntercomparison.hh"
#include "G4INCLClusteringModelNone.hh"

namespace G4INCL {

  namespace Clustering {

    namespace {
      G4ThreadLocal IClusteringModel *theClusteringModel = NULL;
    }

    Cluster* getCluster(Nucleus *n, Particle *p) {
#if !defined(NDEBUG) && !defined(INCLXX_IN_GEANT4_MODE)
      Cluster * const c=theClusteringModel->getCluster(n,p);
// assert(!c || c->getA()<=n->getA()/2);
      return c;
#else
      return theClusteringModel->getCluster(n,p);
#endif
    }

    G4bool clusterCanEscape(Nucleus const * const n, Cluster const * const c) {
      return theClusteringModel->clusterCanEscape(n, c);
    }

    IClusteringModel *getClusteringModel() { return theClusteringModel; }

    void setClusteringModel(IClusteringModel * const model) {
      theClusteringModel = model;
    }

    void deleteClusteringModel() {
      delete theClusteringModel;
      theClusteringModel = NULL;
    }

    void initialize(Config const * const theConfig) {
      ClusterAlgorithmType clusterAlgorithm = theConfig->getClusterAlgorithm();
      if(clusterAlgorithm == IntercomparisonClusterAlgorithm)
        setClusteringModel(new ClusteringModelIntercomparison(theConfig));
      else // if(clusterAlgorithm == NoClusterAlgorithm)
        setClusteringModel(new ClusteringModelNone);
    }

  }
}

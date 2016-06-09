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
 * G4INCLReflectionAvatar.cc
 *
 *  Created on: Jun 8, 2009
 *      Author: Pekka Kaitaniemi
 */

#include "G4INCLSurfaceAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLReflectionChannel.hh"
#include "G4INCLTransmissionChannel.hh"
#include "G4INCLClustering.hh"
#include <sstream>
#include <string>

namespace G4INCL {

  SurfaceAvatar::SurfaceAvatar(G4INCL::Particle *aParticle, G4double time, G4INCL::Nucleus *n)
    :IAvatar(time), theParticle(aParticle), theNucleus(n)
  {
    setType(SurfaceAvatarType);
    // TODO Auto-generated constructor stub
  }

  SurfaceAvatar::~SurfaceAvatar() {

  }

  G4INCL::IChannel* SurfaceAvatar::getChannel() const
  {
    if(!theParticle->isParticipant()) {
      return new ReflectionChannel(theNucleus, theParticle);
    }

    // Don't try to make a cluster if the leading particle is too slow
    const G4double transmissionProbability = theNucleus->getTransmissionProbability(theParticle);

    if(theParticle->isNucleon() && transmissionProbability>1.E-4) {
      Cluster *candidateCluster = 0;

      candidateCluster = Clustering::getCluster(theNucleus, theParticle);
      if(candidateCluster != 0 &&
          Clustering::clusterCanEscape(candidateCluster)) {
        return new TransmissionChannel(theNucleus, candidateCluster);
      } else {
        delete candidateCluster;
      }
    }

    // If we haven't transmitted a cluster (maybe cluster feature was
    // disabled or maybe we just can't produce an acceptable cluster):
    const G4double x = Random::shoot();

    if(x <= transmissionProbability) { // Transmission
      return new TransmissionChannel(theNucleus, theParticle);
    } else { // Reflection
      return new ReflectionChannel(theNucleus, theParticle);
    }
  }

  G4INCL::FinalState* SurfaceAvatar::getFinalState() const
  {
    return getChannel()->getFinalState();
  }

  void SurfaceAvatar::preInteraction() {}
  FinalState *SurfaceAvatar::postInteraction(FinalState *fs) {
    ParticleList outgoing = fs->getOutgoingParticles();
    if(!outgoing.empty()) { // Transmission
      // assert(outgoing.size()==1);
      Particle *out = outgoing.front();
      if(out->isCluster()) {
	Cluster *clusterOut = dynamic_cast<Cluster*>(out);
        ParticleList const *components = clusterOut->getParticles();
        for(ParticleIter i=components->begin(); i!=components->end(); ++i) {
          if((*i)->isParticipant())
            theNucleus->getStore()->getBook()->decrementParticipants();
        }
      } else if(theParticle->isParticipant()) {
        // assert(out==theParticle);
        theNucleus->getStore()->getBook()->decrementParticipants();
      }
    }
    return fs;
  }

  std::string SurfaceAvatar::dump() const {
    std::stringstream ss;
    ss << "(avatar " << theTime << " 'reflection" << std::endl
      << "(list " << std::endl 
      << theParticle->dump()
      << "))" << std::endl;
    return ss.str();
  }
}

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

/*
 * G4INCLBinaryCollisionAvatar.hh
 *
 *  \date Jun 5, 2009
 * \author Pekka Kaitaniemi
 */

#ifndef G4INCLBINARYCOLLISIONAVATAR_HH_
#define G4INCLBINARYCOLLISIONAVATAR_HH_

#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLInteractionAvatar.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {

  class BinaryCollisionAvatar: public InteractionAvatar {
  public:
    BinaryCollisionAvatar(G4double, G4double, G4INCL::Nucleus*, G4INCL::Particle*, G4INCL::Particle*);
    virtual ~BinaryCollisionAvatar();
    G4INCL::IChannel* getChannel();
    ParticleList getParticles() const {
      ParticleList theParticleList;
      theParticleList.push_back(particle1);
      theParticleList.push_back(particle2);
      return theParticleList;
    };

    virtual void preInteraction();
    virtual void postInteraction(FinalState *);

    std::string dump() const;

    static void setCutNN(const G4double c) {
      cutNN = c;
      cutNNSquared = cutNN*cutNN;
    }

    static G4double getCutNN() { return cutNN; }

    static G4double getCutNNSquared() { return cutNNSquared; }
    
    /// \brief Get the bias
    static G4double getBias() { return bias; }
    
    /// \brief Set the bias
    static void setBias(const G4double b) { bias=b; }

  private:
    static G4ThreadLocal G4double cutNN;
    static G4ThreadLocal G4double cutNNSquared;
    static G4ThreadLocal G4double bias;
    G4double theCrossSection;
    G4bool isParticle1Spectator;
    G4bool isParticle2Spectator;
    G4bool isElastic;

    INCL_DECLARE_ALLOCATION_POOL(BinaryCollisionAvatar)
  };

}

#endif /* G4INCLBINARYCOLLISIONAVATAR_HH_ */

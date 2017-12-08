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

#ifndef G4INCLDECAYAVATAR_HH_
#define G4INCLDECAYAVATAR_HH_

#include "G4INCLInteractionAvatar.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {

  /**
   * Decay avatar
   *
   * The reflection avatar is created when a particle reaches the boundary of the nucleus.
   * At this point it can either be reflected from the boundary or exit the nucleus.
   */
  class DecayAvatar: public InteractionAvatar {
  public:
    DecayAvatar(G4INCL::Particle *aParticle, G4double time, G4INCL::Nucleus *aNucleus, G4bool force=false);
    DecayAvatar(G4INCL::Particle *aParticle, G4INCL::Particle *bParticle, G4double time, G4INCL::Nucleus *aNucleus, G4bool force=false);
    virtual ~DecayAvatar();

    IChannel* getChannel();
    void fillFinalState(FinalState *fs);

    virtual void preInteraction();
    virtual void postInteraction(FinalState *fs);

    ParticleList getParticles() const {
      ParticleList theParticleList;
      theParticleList.push_back(particle1);
      return theParticleList;
    }

    std::string dump() const;
  private:
    G4bool forced;
    ThreeVector const incidentDirection;

    INCL_DECLARE_ALLOCATION_POOL(DecayAvatar)
  };

}

#endif /* G4INCLDECAYAVATAR_HH_ */

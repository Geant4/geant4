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

/*
 * G4INCLSurfaceAvatar.hh
 *
 *  \date Jun 8, 2009
 * \author Pekka Kaitaniemi
 */

#ifndef G4INCLSURFACEAVATAR_HH_
#define G4INCLSURFACEAVATAR_HH_

#include "G4INCLIAvatar.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"

namespace G4INCL {

  /**
   * Surface avatar
   *
   * The reflection avatar is created when a particle reaches the boundary of the nucleus.
   * At this point it can either be reflected from the boundary or exit the nucleus.
   */
  class SurfaceAvatar: public G4INCL::IAvatar {
  public:
    SurfaceAvatar(G4INCL::Particle *aParticle, G4double time, G4INCL::Nucleus *aNucleus);
    virtual ~SurfaceAvatar();

    G4INCL::IChannel* getChannel() const;
    G4INCL::FinalState* getFinalState() const;

    virtual void preInteraction();
    virtual FinalState *postInteraction(FinalState *);

    ParticleList getParticles() const {
      ParticleList theParticleList;
      theParticleList.push_back(theParticle);
      return theParticleList;
    }

    std::string dump() const;

    /// \brief Calculate the transmission probability for the particle
    G4double getTransmissionProbability(Particle const * const particle) const;

  private:
    G4INCL::Particle *theParticle;
    G4INCL::Nucleus *theNucleus;
  };

}

#endif /* G4INCLSURFACEAVATAR_HH_ */

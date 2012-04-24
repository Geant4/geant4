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
// INCL++ revision: v5.0.5
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4INCLParticleEntryChannel.hh"
#include "G4INCLRootFinder.hh"

namespace G4INCL {

  ParticleEntryChannel::ParticleEntryChannel(Nucleus *n, Particle *p)
    :theNucleus(n), theParticle(p)
  {}

  ParticleEntryChannel::~ParticleEntryChannel()
  {}

  FinalState* ParticleEntryChannel::getFinalState() {
    const G4double energyBefore = theParticle->getEnergy();
    particleEnters();
    theNucleus->insertParticipant(theParticle);

    FinalState *fs = new FinalState();
    fs->addModifiedParticle(theParticle);
    fs->setTotalEnergyBeforeInteraction(energyBefore);
    return fs;
  }

  void ParticleEntryChannel::particleEnters() {

    // TODO: this is the place to add refraction

    // Add the nuclear potential to the kinetic energy when entering the
    // nucleus

    class IncomingEFunctor : public RootFunctor {
      public:
        IncomingEFunctor(Particle * const p, NuclearPotential::INuclearPotential const * const np) :
          theParticle(p), thePotential(np) {
            theEnergy=theParticle->getEnergy();
          }
        ~IncomingEFunctor() {}
        G4double operator()(const G4double v) const {
          theParticle->setEnergy(theEnergy + v);
          theParticle->setPotentialEnergy(v);
          // Scale the particle momentum
          theParticle->adjustMomentumFromEnergy();
          return v - thePotential->computePotentialEnergy(theParticle);
        }
        void cleanUp(const G4bool /*success*/) const {}
      private:
        Particle *theParticle;
        NuclearPotential::INuclearPotential const *thePotential;
        G4double theEnergy;
    } theIncomingEFunctor(theParticle,theNucleus->getPotential());

    G4double v = theNucleus->getPotential()->computePotentialEnergy(theParticle);
    G4bool success = RootFinder::solve(&theIncomingEFunctor, v);
    if(success) { // Apply the solution
      std::pair<G4double,G4double> theSolution = RootFinder::getSolution();
      theIncomingEFunctor(theSolution.first);
    } else {
      WARN("Couldn't compute the potential for incoming particle, root-finding algorithm failed." << std::endl);
    }
  }

}


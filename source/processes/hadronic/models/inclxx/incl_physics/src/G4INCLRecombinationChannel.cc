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

/** \file G4INCLRecombinationChannel.cc
 * \brief Delta-nucleon recombination channel.
 *
 *  \date 25 March 2011
 * \author Davide Mancusi
 */

#include "G4INCLRecombinationChannel.hh"
#include "G4INCLRandom.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLGlobals.hh"

// #include <cassert>

namespace G4INCL {

  RecombinationChannel::RecombinationChannel(Particle *p1, Particle *p2)
  {
    if(p1->isDelta()) {
// assert(p2->isNucleon());
      theDelta = p1;
      theNucleon = p2;
    } else {
// assert(p1->isNucleon());
      theDelta = p2;
      theNucleon = p1;
    }
  }

  RecombinationChannel::~RecombinationChannel()
  {
  }

  void RecombinationChannel::fillFinalState(FinalState *fs)
  {
    // Compute the total available energy in the CM
    const G4double sqrts = KinematicsUtils::totalEnergyInCM(theDelta, theNucleon);

    // Assign the types of the final-state particles
    switch(theDelta->getType()) {
      case DeltaPlusPlus:
// assert(theNucleon->getType()!=Proton);
        theDelta->setType(Proton);
        theNucleon->setType(Proton);
        break;
      case DeltaPlus:
        theDelta->setType(Proton);
        break;
      case DeltaZero:
        theDelta->setType(Neutron);
        break;
      case DeltaMinus:
// assert(theNucleon->getType()!=Neutron);
        theDelta->setType(Neutron);
        theNucleon->setType(Neutron);
        break;
      default:
        INCL_ERROR("Unknown particle type in RecombinationChannel" << '\n');
        break;
    }

    // Calculate the momenta of the nucleons in the final state
    const G4double pCM = KinematicsUtils::momentumInCM(sqrts, theDelta->getMass(), theNucleon->getMass());

    // The angular distribution of final-state nucleons is isotropic
    ThreeVector momentum = Random::normVector(pCM);

    // Assign the momenta
    theDelta->setMomentum(momentum);
    theNucleon->setMomentum(-momentum);

    // Update the kinetic energies
    theDelta->adjustEnergyFromMomentum();
    theNucleon->adjustEnergyFromMomentum();

    // Create the final state
    fs->addModifiedParticle(theDelta);
    fs->addModifiedParticle(theNucleon);

  }

}

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

#include "G4INCLEtaNToPiPiNChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {

  const G4double EtaNToPiPiNChannel::angularSlope = 15.; // from piNToMultipions

  EtaNToPiPiNChannel::EtaNToPiPiNChannel(Particle *p1, Particle *p2)
    : ind2(0),
    particle1(p1),
    particle2(p2)
  {}

  EtaNToPiPiNChannel::~EtaNToPiPiNChannel(){

  }

  void EtaNToPiPiNChannel::fillFinalState(FinalState *fs) {

    Particle * nucleon;
    Particle * eta;
    if(particle1->isNucleon()) {
      nucleon = particle1;
      eta = particle2;
    } else {
      nucleon = particle2;
      eta = particle1;
    }

    const G4double sqrtS = KinematicsUtils::totalEnergyInCM(nucleon, eta);
			
    ind2=ParticleTable::getIsospin(nucleon->getType());

    eta->setType(PiZero);

    ParticleType pionType;
    G4double rs1=Random::shoot();
    if (ind2 == 1) {
        if (rs1*6. > 5.) pionType=PiZero;
        else if (rs1*6. > 3.) {
            pionType=PiPlus;
            ind2=-ind2;
        }
        else {
            pionType=PiPlus;
            eta->setType(PiMinus);
            }
    }
    else {
        if (rs1*6. > 5.) pionType=PiZero;
        else if (rs1*6. > 3.) {
            pionType=PiMinus;
            ind2=-ind2;
        }
        else {
            pionType=PiPlus;
            eta->setType(PiMinus);
        }
    }
			

    const ParticleType tn=ParticleTable::getNucleonType(ind2);
    nucleon->setType(tn);
    ParticleList list;
    list.push_back(nucleon);
    list.push_back(eta);
    const ThreeVector &rcolpion = eta->getPosition();
    const ThreeVector zero;
    Particle *newPion = new Particle(pionType,zero,rcolpion);
    list.push_back(newPion);
    fs->addModifiedParticle(nucleon);
    fs->addModifiedParticle(eta);
    fs->addCreatedParticle(newPion);

    PhaseSpaceGenerator::generateBiased(sqrtS, list, 0, angularSlope);

  }

}

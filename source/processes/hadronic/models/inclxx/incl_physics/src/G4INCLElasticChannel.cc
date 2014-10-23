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

#include "G4INCLElasticChannel.hh"
#include "G4INCLRandom.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLCrossSections.hh"
#include "G4INCLGlobals.hh"

namespace G4INCL {

  ElasticChannel::ElasticChannel(Particle *p1, Particle *p2)
    :particle1(p1), particle2(p2)
  {
  }

  ElasticChannel::~ElasticChannel()
  {
  }

  void ElasticChannel::fillFinalState(FinalState *fs)
  {
    ParticleType p1TypeOld = particle1->getType();
    ParticleType p2TypeOld = particle2->getType();

    /* Concerning the way we calculate the lab momentum, see the considerations
     * in CrossSections::elasticNNLegacy().
     */
    const G4double s = KinematicsUtils::squareTotalEnergyInCM(particle1, particle2);
    const G4double pl = KinematicsUtils::momentumInLab(s, ParticleTable::effectiveNucleonMass, ParticleTable::effectiveNucleonMass);

    const G4int isospin = ParticleTable::getIsospin(particle1->getType()) +
      ParticleTable::getIsospin(particle2->getType());

    // Calculate the outcome of the channel:
    G4double psq = particle1->getMomentum().mag2();
    G4double pnorm = std::sqrt(psq);
    G4double b = CrossSections::calculateNNAngularSlope(pl, isospin);
    G4double btmax = 4.0 * psq * b;
    G4double z = std::exp(-btmax);
    G4double ranres = Random::shoot();
    G4double y = 1.0 - ranres * (1.0 - z);
    G4double T = std::log(y)/b;
    G4int iexpi = 0;
    G4double apt = 1.0;

    // Handle np case
    if((particle1->getType() == Proton && particle2->getType() == Neutron) ||
       (particle1->getType() == Neutron && particle2->getType() == Proton)) {
      if(pl > 800.0) {
        const G4double x = 0.001 * pl; // Transform to GeV
        apt = (800.0/pl)*(800.0/pl);
        G4double cpt = std::max(6.23 * std::exp(-1.79*x), 0.3);
        G4double alphac = 100.0 * 1.0e-6;
        G4double aaa = (1 + apt) * (1 - std::exp(-btmax))/b;
        G4double argu = psq * alphac;

        if(argu >= 8) {
          argu = 0.0;
        } else {
          argu = std::exp(-4.0 * argu);
        }

        G4double aac = cpt * (1.0 - argu)/alphac;
        G4double fracpn = aaa/(aac + aaa);
        if(Random::shoot() > fracpn) {
          z = std::exp(-4.0 * psq *alphac);
          iexpi = 1;
          y = 1.0 - ranres*(1.0 - z);
          T = std::log(y)/alphac;
        }
      }
    }

    G4double ctet = 1.0 + 0.5*T/psq;
    if(std::abs(ctet) > 1.0) ctet = Math::sign(ctet);
    G4double stet = std::sqrt(1.0 - ctet*ctet);
    G4double rndm = Random::shoot();

    G4double fi = Math::twoPi * rndm;
    G4double cfi = std::cos(fi);
    G4double sfi = std::sin(fi);

    G4double xx = particle1->getMomentum().perp2();
    G4double zz = std::pow(particle1->getMomentum().getZ(), 2);

    if(xx >= (zz * 1.0e-8)) {
      ThreeVector p = particle1->getMomentum();
      G4double yn = std::sqrt(xx);
      G4double zn = yn * pnorm;
      G4double ex[3], ey[3], ez[3];
      ez[0] = p.getX() / pnorm;
      ez[1] = p.getY() / pnorm;
      ez[2] = p.getZ() / pnorm;

      // Vector Ex is chosen arbitrarily:
      ex[0] = p.getY() / yn;
      ex[1] = -p.getX() / yn;
      ex[2] = 0.0;

      ey[0] = p.getX() * p.getZ() / zn;
      ey[1] = p.getY() * p.getZ() / zn;
      ey[2] = -xx/zn;

      G4double pX = (ex[0]*cfi*stet + ey[0]*sfi*stet + ez[0]*ctet) * pnorm;
      G4double pY = (ex[1]*cfi*stet + ey[1]*sfi*stet + ez[1]*ctet) * pnorm;
      G4double pZ = (ex[2]*cfi*stet + ey[2]*sfi*stet + ez[2]*ctet) * pnorm;

      ThreeVector p1momentum = ThreeVector(pX, pY, pZ);
      particle1->setMomentum(p1momentum);
      particle2->setMomentum(-p1momentum);
    } else { // if(xx < (zz * 1.0e-8)) {
      G4double momZ = particle1->getMomentum().getZ();
      G4double pX = momZ * cfi * stet;
      G4double pY = momZ * sfi * stet;
      G4double pZ = momZ * ctet;

      ThreeVector p1momentum(pX, pY, pZ);
      particle1->setMomentum(p1momentum);
      particle2->setMomentum(-p1momentum);
    }

    // Handle backward scattering here.

    if((particle1->getType() == Proton && particle2->getType() == Neutron) ||
       (particle1->getType() == Neutron && particle2->getType() == Proton)) {
      rndm = Random::shoot();
      apt = 1.0;
      if(pl > 800.0) {
        apt = std::pow(800.0/pl, 2);
      }
      if(iexpi == 1 || rndm > 1.0/(1.0 + apt)) {
        particle1->setType(p2TypeOld);
        particle2->setType(p1TypeOld);
      }
    }

    // Note: there is no need to update the kinetic energies of the particles,
    // as this is elastic scattering.

    fs->addModifiedParticle(particle1);
    fs->addModifiedParticle(particle2);

    }

}

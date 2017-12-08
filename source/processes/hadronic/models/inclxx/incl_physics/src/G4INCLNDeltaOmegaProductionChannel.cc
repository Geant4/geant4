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

#include "G4INCLNDeltaOmegaProductionChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {

  const G4double NDeltaOmegaProductionChannel::angularSlope = 6.;
  const G4int NDeltaOmegaProductionChannel::maxTries = 100000;

  NDeltaOmegaProductionChannel::NDeltaOmegaProductionChannel(Particle *p1,Particle *p2)
    : particle1(p1), particle2(p2)
  {}

  NDeltaOmegaProductionChannel::~NDeltaOmegaProductionChannel() {}

  G4double NDeltaOmegaProductionChannel::sampleDeltaMass(G4double ecmorigin) {
    const G4double ecm = ecmorigin - 783.437; // 783.437 MeV translation to open pion(delta) production in NNOmega
    const G4double maxDeltaMass = ecm - ParticleTable::effectiveNucleonMass - 1.0;
    const G4double maxDeltaMassRndm = std::atan((maxDeltaMass-ParticleTable::effectiveDeltaMass)*2./ParticleTable::effectiveDeltaWidth);
    const G4double deltaMassRndmRange = maxDeltaMassRndm - ParticleTable::minDeltaMassRndm;
// assert(deltaMassRndmRange>0.);

    G4double y=ecm*ecm;
    G4double q2=(y-1.157776E6)*(y-6.4E5)/y/4.0; // 1.157776E6 = 1076^2, 6.4E5 = 800^2
    G4double q3=std::pow(std::sqrt(q2), 3.);
    const G4double f3max=q3/(q3+5.832E6); // 5.832E6 = 180^3
    G4double x;

    G4int nTries = 0;
    G4bool success = false;
    while(!success) { /* Loop checking, 10.07.2015, D.Mancusi */
      if(++nTries >= maxTries) {
        INCL_WARN("NDeltaOmegaProductionChannel::sampleDeltaMass loop was stopped because maximum number of tries was reached. Minimum delta mass "
                  << ParticleTable::minDeltaMass << " MeV with CM energy " << ecm << " MeV may be unphysical." << '\n');
        return ParticleTable::minDeltaMass;
      }

      G4double rndm = ParticleTable::minDeltaMassRndm + Random::shoot() * deltaMassRndmRange;
      y = std::tan(rndm);
      x = ParticleTable::effectiveDeltaMass + 0.5*ParticleTable::effectiveDeltaWidth*y;
// assert(x>=ParticleTable::minDeltaMass && ecm >= x + ParticleTable::effectiveNucleonMass + 1.0);

      // generation of the delta mass with the penetration factor
      // (see prc56(1997)2431)
      y=x*x;
      q2=(y-1.157776E6)*(y-6.4E5)/y/4.0; // 1.157776E6 = 1076^2, 6.4E5 = 800^2
      q3=std::pow(std::sqrt(q2), 3.);
      const G4double f3=q3/(q3+5.832E6); // 5.832E6 = 180^3
      rndm = Random::shoot();
      if (rndm*f3max < f3)
        success = true;
    }
    return x;
  }

 void NDeltaOmegaProductionChannel::fillFinalState(FinalState *fs) {
  
/**
*
* Unlike NN -> NDelta, NN -> NDeltaOmega is drawn from a phase-space generator
*
**/

  G4int is1=ParticleTable::getIsospin(particle1->getType());
  G4int is2=ParticleTable::getIsospin(particle2->getType());
  
  ParticleList list;
  list.push_back(particle1);
  list.push_back(particle2);
  
//  isospin Repartition of N and Delta;
  G4double ecm = KinematicsUtils::totalEnergyInCM(particle1, particle2);  
  const G4int isospin = is1+is2;
  
  G4double rndm = 0.0;
  G4double xmdel = sampleDeltaMass(ecm);

  G4int index2=0;
  if (isospin == 0) { // pn case
   rndm = Random::shoot();
   if (rndm < 0.5) index2=1;
  }

  if (isospin == 0) {
   if(index2 == 1) {
    G4int isi=is1;
    is1=is2;
    is2=isi;
   }
//   particle1->setHelicity(0.0);
  } else {
   rndm = Random::shoot();
   if (rndm >= 0.25) {
    is1=3*is1;
    is2=-is2;
   }
//   particle1->setHelicity(ctet*ctet);
  }
  
  if(is1 == ParticleTable::getIsospin(DeltaMinus)) {
   particle1->setType(DeltaMinus);
  } else if(is1 == ParticleTable::getIsospin(DeltaZero)) {
   particle1->setType(DeltaZero);
  } else if(is1 == ParticleTable::getIsospin(DeltaPlus)) {
   particle1->setType(DeltaPlus);
  } else if(is1 == ParticleTable::getIsospin(DeltaPlusPlus)) {
   particle1->setType(DeltaPlusPlus);
  }
  
  if(is2 == ParticleTable::getIsospin(Proton)) {
   particle2->setType(Proton);
  } else if(is2 == ParticleTable::getIsospin(Neutron)) {
   particle2->setType(Neutron);
  }
  
  if(particle1->isDelta()) particle1->setMass(xmdel);
  if(particle2->isDelta()) particle2->setMass(xmdel);
  
  
  const ThreeVector &rcolnucleon1 = particle1->getPosition();
  const ThreeVector &rcolnucleon2 = particle2->getPosition();
  const ThreeVector rcol = (rcolnucleon1+rcolnucleon2)*0.5;
  const ThreeVector zero;
  Particle *omega = new Particle(Omega,zero,rcol);
  list.push_back(omega);
  fs->addCreatedParticle(omega);
  
  const G4double sqrtS = KinematicsUtils::totalEnergyInCM(particle1, particle2);
  G4int biasIndex = ((Random::shoot()<0.5) ? 0 : 1);
  PhaseSpaceGenerator::generateBiased(sqrtS, list, biasIndex, angularSlope);
  
  const ThreeVector vz(0.0,0.0,1.0);
  G4double ctet=(particle1->getMomentum().dot(vz))/particle1->getMomentum().mag();
  if (isospin == 0)
   particle1->setHelicity(0.0);
  else
   particle1->setHelicity(ctet*ctet);
  fs->addModifiedParticle(particle1);
  fs->addModifiedParticle(particle2);
  
 }
 
}

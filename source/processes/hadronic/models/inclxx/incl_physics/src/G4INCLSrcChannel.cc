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
 * \file G4INCLSrcChannel.cc
 *
 * \date Feb 24, 2022
 * \author Jose Luis Rodriguez-Sanchez
 */

#include "G4INCLSrcChannel.hh"
#include "G4INCLCrossSections.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLInteractionAvatar.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLRandom.hh"

namespace G4INCL {

SrcChannel::SrcChannel(Particle *p1, Particle *p2, Nucleus *n)
    : particle1(p1), particle2(p2), thenucleus(n) {
  fDistSrc = ParticleTable::getsrcPairDistance();
  srcpartner = nullptr;
  ftype1 = UnknownParticle;
  ftype2 = UnknownParticle;
}

SrcChannel::~SrcChannel() {}

Particle *SrcChannel::findpairpartner(Particle *pt) {

  const auto pair = pt->getSrcPair();
  const ParticleType t = pt->getType();

  ParticleList const &inside = thenucleus->getStore()->getParticles();
  for (ParticleIter p = inside.begin(), e = inside.end(); p != e; ++p) {
    if ((*p)->getSrcPair() == pair &&
        (t != (*p)->getType() ||
         (t == (*p)->getType() && pt->getID() != (*p)->getID()))) 
         {
          return (*p);
         }
  }
  INCL_ERROR("SrcChannel: pair not found" << '\n');
  return NULL;
}

void SrcChannel::fillFinalState(FinalState *fs, ParticleType type1,
                                ParticleType type2) {
  ftype1 = type1;
  ftype2 = type2;
  fillFinalState(fs);
}

void SrcChannel::fillFinalState(FinalState *fs) {

  auto fSource = InteractionAvatar::Instance();

  G4double psrcmax = 2.0 * PhysicalConstants::Pf; // Fermi momentum in MeV/c

  if (particle1->getSrcPair() > 0) {

    if (particle2->getSrcPair() > 0) 
    {
      INCL_ERROR(particle2->print() << " \n");
    }
    srcpartner = findpairpartner(particle1);

    if (srcpartner) 
    {
      srcpartner->setSrcPartner();
      fSource->setSrcPartner(srcpartner);

      ThreeVector d1 = particle1->getPosition();
      ThreeVector d2 = srcpartner->getPosition();
      auto d = d1 - d2;
      if (d.mag() > fDistSrc) {
        INCL_DEBUG("Distance src > " << fDistSrc << " fm : " << d.mag()
                                     << " \n");
      }

      auto x = (fDistSrc - d.mag()) / fDistSrc;
      auto srcp = x * x * psrcmax;
      INCL_DEBUG("Src momentum = " << srcp << " , eventnb: "
                                   << theEventInfo.eventNumber << " \n");

      auto pmomentum =
          particle1->getMomentum() / particle1->getMomentum().mag() * srcp;
      particle1->setMomentum(particle1->getMomentum() + pmomentum);
      srcpartner->setMomentum(srcpartner->getMomentum() - pmomentum);

      particle1->adjustEnergyFromMomentum();
      srcpartner->adjustEnergyFromMomentum();
      thenucleus->updatePotentialEnergy(srcpartner);

      fs->addModifiedParticle(particle1);
      fs->addModifiedParticle(particle2);
      fs->addModifiedParticle(srcpartner);
    }
  } else {

    if (particle1->getSrcPair() > 0) 
    {
      INCL_ERROR(particle1->print() << " \n");
    }
    srcpartner = findpairpartner(particle2);

    if (srcpartner) 
    {
      srcpartner->setSrcPartner();
      fSource->setSrcPartner(srcpartner);

      ThreeVector d1 = particle2->getPosition();
      ThreeVector d2 = srcpartner->getPosition();
      auto d = d1 - d2;

      if (d.mag() > fDistSrc) {
        INCL_DEBUG("Distance src > " << fDistSrc << " fm : " << d.mag()
                                     << " \n");
      }

      auto x = (fDistSrc - d.mag()) / fDistSrc;
      auto srcp = x * x * psrcmax;
      INCL_DEBUG("Src momentum = " << srcp << " , eventnb: "
                                   << theEventInfo.eventNumber << " \n");

      auto pmomentum =
          particle2->getMomentum() / particle2->getMomentum().mag() * srcp;

      particle2->setMomentum(particle2->getMomentum() + pmomentum);
      srcpartner->setMomentum(srcpartner->getMomentum() - pmomentum);

      particle2->adjustEnergyFromMomentum();
      srcpartner->adjustEnergyFromMomentum();
      thenucleus->updatePotentialEnergy(srcpartner);

      fs->addModifiedParticle(particle1);
      fs->addModifiedParticle(particle2);
      fs->addModifiedParticle(srcpartner);
    }
  }
  thenucleus->getStore()->getBook().incrementAcceptedSrcCollisions();
}
} // namespace G4INCL

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

#include "G4INCLFinalState.hh"

namespace G4INCL {

  FinalState::FinalState() :
    totalEnergyBeforeInteraction(0.0), validity(ValidFS),
    blockedDelta(NULL)
  {
  }

  FinalState::~FinalState()
  {
  }

  void FinalState::addModifiedParticle(Particle *p)
  {
    modified.push_back(p);
  }

  void FinalState::addOutgoingParticle(Particle *p)
  {
    outgoing.push_back(p);
  }

  void FinalState::addDestroyedParticle(Particle *p)
  {
    destroyed.push_back(p);
  }

  void FinalState::addCreatedParticle(Particle *p)
  {
    created.push_back(p);
  }

  void FinalState::addEnteringParticle(Particle *p)
  {
    entering.push_back(p);
  }

  ParticleList const &FinalState::getModifiedParticles() const
  {
    return modified;
  }

  ParticleList const &FinalState::getOutgoingParticles() const
  {
    return outgoing;
  }

  ParticleList const &FinalState::getDestroyedParticles() const
  {
    return destroyed;
  }

  ParticleList const &FinalState::getCreatedParticles() const
  {
    return created;
  }

  ParticleList const &FinalState::getEnteringParticles() const
  {
    return entering;
  }

  std::string FinalState::print() const {
    std::stringstream ss;
    ss << "Modified particles:" << std::endl;
    for(ParticleIter iter = modified.begin(); iter != modified.end(); ++iter)
      ss << (*iter)->print();
    ss << "Outgoing particles:" << std::endl;
    for(ParticleIter iter = outgoing.begin(); iter != outgoing.end(); ++iter)
      ss << (*iter)->print();
    ss << "Destroyed particles:" << std::endl;
    for(ParticleIter iter = destroyed.begin(); iter != destroyed.end(); ++iter)
      ss << (*iter)->print();
    ss << "Created particles:" << std::endl;
    for(ParticleIter iter = created.begin(); iter != created.end(); ++iter)
      ss << (*iter)->print();
    ss << "Entering particles:" << std::endl;
    for(ParticleIter iter = entering.begin(); iter != entering.end(); ++iter)
      ss << (*iter)->print();
    return ss.str();
  }

}

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

#include "G4INCLFinalState.hh"

namespace G4INCL {

  FinalState::FinalState() {
    reset();
  }

  FinalState::~FinalState()
  {
  }

  void FinalState::reset() {
    totalEnergyBeforeInteraction = 0.0;
    validity = ValidFS;
    outgoing.clear();
    created.clear();
    destroyed.clear();
    modified.clear();
    entering.clear();
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
    ss << "Modified particles:" << '\n';
    for(ParticleIter iter=modified.begin(), e=modified.end(); iter!=e; ++iter)
      ss << (*iter)->print();
    ss << "Outgoing particles:" << '\n';
    for(ParticleIter iter=outgoing.begin(), e=outgoing.end(); iter!=e; ++iter)
      ss << (*iter)->print();
    ss << "Destroyed particles:" << '\n';
    for(ParticleIter iter=destroyed.begin(), e=destroyed.end(); iter!=e; ++iter)
      ss << (*iter)->print();
    ss << "Created particles:" << '\n';
    for(ParticleIter iter=created.begin(), e=created.end(); iter!=e; ++iter)
      ss << (*iter)->print();
    ss << "Entering particles:" << '\n';
    for(ParticleIter iter=entering.begin(), e=entering.end(); iter!=e; ++iter)
      ss << (*iter)->print();
    return ss.str();
  }

}

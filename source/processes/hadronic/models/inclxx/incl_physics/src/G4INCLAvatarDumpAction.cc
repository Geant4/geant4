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

#include "G4INCLAvatarDumpAction.hh"
#include <sstream>
#include <string>

namespace G4INCL {

  AvatarDumpAction::AvatarDumpAction() :
    oFile(0),
    eventCounter(0)
  {
  }

  AvatarDumpAction::~AvatarDumpAction() {}

  void AvatarDumpAction::beforeCascadeUserAction(IPropagationModel * /*pm*/) {
    std::stringstream ss;
    ss << "avatar-dump-" << eventCounter << ".dat";
    oFile = new std::ofstream(ss.str().c_str());
  }

  void AvatarDumpAction::afterAvatarUserAction(IAvatar *avatar, Nucleus *nucleus, FinalState *finalState) {
    ParticleList particles = nucleus->getStore()->getParticles();
    ParticleList highlight;
    if(finalState) {
      ParticleList const &modified = finalState->getModifiedParticles();
      highlight.insert(highlight.end(), modified.begin(), modified.end());
      ParticleList const &outgoing = finalState->getOutgoingParticles();
      highlight.insert(highlight.end(), outgoing.begin(), outgoing.end());
      ParticleList const &destroyed = finalState->getDestroyedParticles();
      highlight.insert(highlight.end(), destroyed.begin(), destroyed.end());
      ParticleList const &created = finalState->getCreatedParticles();
      highlight.insert(highlight.end(), created.begin(), created.end());
      ParticleList const &entering = finalState->getEnteringParticles();
      highlight.insert(highlight.end(), entering.begin(), entering.end());
      particles.insert(particles.end(), created.begin(), created.end());
      particles.insert(particles.end(), entering.begin(), entering.end());
    }

    (*oFile) << avatar->getTime() << '\t' << avatar->getType() << '\t' << particles.size() << '\n';
    for(ParticleIter p=particles.begin(), e=particles.end(); p!=e; ++p) {
      ThreeVector const &pos = (*p)->getPosition();
      ThreeVector const &vel = (*p)->getPropagationVelocity();
      G4int highlightIt = highlight.contains(*p);
      (*oFile)
        << (*p)->getID() << '\t'
        << (*p)->getParticipantType() << '\t'
        << (*p)->getType() << '\t'
        << pos.getX() << '\t'
        << pos.getY() << '\t'
        << pos.getZ() << '\t'
        << vel.getX() << '\t'
        << vel.getY() << '\t'
        << vel.getZ() << '\t'
        << (*p)->getKineticEnergy() << '\t'
        << (*p)->getPotentialEnergy() << '\t'
        << highlightIt << '\n';
    }
  }

  void AvatarDumpAction::afterCascadeUserAction(Nucleus * /*nucleus*/) {
    oFile->close();
    delete oFile;
    ++eventCounter;
  }

}

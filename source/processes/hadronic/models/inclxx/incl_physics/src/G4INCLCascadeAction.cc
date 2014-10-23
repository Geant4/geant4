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

#include "G4INCLCascadeAction.hh"
#include "G4INCLLogger.hh"
#include "G4INCLRandom.hh"

namespace G4INCL {

  CascadeAction::CascadeAction() :
    stepCounter(0)
  {}

  CascadeAction::~CascadeAction()
  {}

  void CascadeAction::beforeRunAction(Config const *config) {
    beforeRunDefaultAction(config);
    beforeRunUserAction(config);
  }

  void CascadeAction::beforeCascadeAction(IPropagationModel *pm) {
    beforeCascadeDefaultAction(pm);
    beforeCascadeUserAction(pm);
  }

  void CascadeAction::beforePropagationAction(IPropagationModel *pm) {
    beforePropagationDefaultAction(pm);
    beforePropagationUserAction(pm);
  }

  void CascadeAction::beforeAvatarAction(IAvatar *a, Nucleus *n) {
    beforeAvatarDefaultAction(a, n);
    beforeAvatarUserAction(a, n);
  }

  void CascadeAction::afterAvatarAction(IAvatar *a, Nucleus *n, FinalState *fs) {
    afterAvatarDefaultAction(a, n, fs);
    afterAvatarUserAction(a, n, fs);
  }

  void CascadeAction::afterPropagationAction(IPropagationModel *pm, IAvatar *avatar) {
    afterPropagationDefaultAction(pm, avatar);
    afterPropagationUserAction(pm, avatar);
  }

  void CascadeAction::afterCascadeAction(Nucleus *n) {
    afterCascadeDefaultAction(n);
    afterCascadeUserAction(n);
  }

  void CascadeAction::afterRunAction() {
    afterRunDefaultAction();
    afterRunUserAction();
  }



  void CascadeAction::beforeRunDefaultAction(Config const * /*pm*/) {}

  void CascadeAction::beforeCascadeDefaultAction(IPropagationModel * /*pm*/) {}

  void CascadeAction::beforePropagationDefaultAction(IPropagationModel * /*pm*/) {
    // assert(pm->getNucleus()->getStore()->getBook().getCascading() == pm->getNucleus()->getStore()->countCascading());
  }

  void CascadeAction::beforeAvatarDefaultAction(IAvatar *a, Nucleus *n) {
    n->getStore()->getBook().incrementAvatars(a->getType());
    INCL_DEBUG("Random seeds before avatar " << a->getID() << ": "
          << G4INCL::Random::getSeeds() << '\n');
    INCL_DEBUG("Next avatar:" << '\n' << a->dump() << '\n');
  }

  void CascadeAction::afterAvatarDefaultAction(IAvatar *a, Nucleus * /*n*/, FinalState *fs) {

    if(!fs) // do nothing if there is no final state
      return;

    INCL_DEBUG("Random seeds after avatar " << a->getID() << ": "
          << G4INCL::Random::getSeeds() << '\n');

    ParticleList const &modified = fs->getModifiedParticles();
    for(ParticleIter p=modified.begin(), e=modified.end(); p!=e; ++p )
      if(a->isACollision())
        (*p)->incrementNumberOfCollisions();
      else if(a->isADecay())
        (*p)->incrementNumberOfDecays();

    ParticleList const &created = fs->getCreatedParticles();
    for(ParticleIter p=created.begin(), e=created.end(); p!=e; ++p )
      if(a->isACollision())
        (*p)->incrementNumberOfCollisions();
      else if(a->isADecay())
        (*p)->incrementNumberOfDecays();

  }

  void CascadeAction::afterPropagationDefaultAction(IPropagationModel * /* pm */,
                                                IAvatar * /*avatar */) {
    ++stepCounter; // Increment the step counter

#ifdef INCL_DEBUG_LOG
    //   INCL_DATABLOCK(pm->getNucleus()->getStore()->printParticleConfiguration());
#endif
  }

  void CascadeAction::afterCascadeDefaultAction(Nucleus * /*pm*/) {}

  void CascadeAction::afterRunDefaultAction() {}

}

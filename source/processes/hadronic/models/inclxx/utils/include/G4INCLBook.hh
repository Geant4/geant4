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

#ifndef G4INCLBook_hh
#define G4INCLBook_hh 1

#include <map>
#include "G4INCLIAvatar.hh"

namespace G4INCL {
  class Book {
  public:
    Book() {
      reset();
    }
    ~Book() {};

    void reset() {
      nAcceptedCollisions = 0;
      nBlockedCollisions = 0;
      nAcceptedDecays = 0;
      nBlockedDecays = 0;
      currentTime = 0.0;
      firstCollisionTime = 0.0;
      firstCollisionXSec = 0.0;
      firstCollisionSpectatorPosition = 0.0;
      firstCollisionSpectatorMomentum = 0.0;
      firstCollisionIsElastic = false;
      nAvatars[SurfaceAvatarType] = 0;
      nAvatars[CollisionAvatarType] = 0;
      nAvatars[DecayAvatarType] = 0;
      nAvatars[ParticleEntryAvatarType] = 0;
      nCascading = 0;
      nEmittedClusters = 0;
      nEnergyViolationInteraction = 0;
    };

    void incrementAcceptedCollisions() { nAcceptedCollisions++; };
    void incrementBlockedCollisions() { nBlockedCollisions++; };
    void incrementAcceptedDecays() { nAcceptedDecays++; };
    void incrementBlockedDecays() { nBlockedDecays++; };
    void incrementAvatars(AvatarType type) { nAvatars[type]++; };
    void incrementCascading() { nCascading++; }
    void decrementCascading() { nCascading--; }
    void incrementEmittedClusters() { nEmittedClusters++; }
    void incrementEnergyViolationInteraction() { nEnergyViolationInteraction++; }

    void setFirstCollisionTime(const G4double t) { firstCollisionTime = t; };
    G4double getFirstCollisionTime() const { return firstCollisionTime; };

    void setFirstCollisionXSec(const G4double x) { firstCollisionXSec = x; };
    G4double getFirstCollisionXSec() const { return firstCollisionXSec; };

    void setFirstCollisionSpectatorPosition(const G4double x) { firstCollisionSpectatorPosition = x; };
    G4double getFirstCollisionSpectatorPosition() const { return firstCollisionSpectatorPosition; };

    void setFirstCollisionSpectatorMomentum(const G4double x) { firstCollisionSpectatorMomentum = x; };
    G4double getFirstCollisionSpectatorMomentum() const { return firstCollisionSpectatorMomentum; };

    void setFirstCollisionIsElastic(const G4bool e) { firstCollisionIsElastic = e; };
    G4bool getFirstCollisionIsElastic() const { return firstCollisionIsElastic; };

    void setCurrentTime(G4double t) { currentTime = t; };
    G4double getCurrentTime() const { return currentTime; };

    G4int getAcceptedCollisions() const { return nAcceptedCollisions; };
    G4int getBlockedCollisions() const {return nBlockedCollisions; };
    G4int getAcceptedDecays() const { return nAcceptedDecays; };
    G4int getBlockedDecays() const {return nBlockedDecays; };
    G4int getAvatars(AvatarType type) const { return nAvatars.find(type)->second; };
    G4int getCascading() const { return nCascading; };
    G4int getEmittedClusters() const { return nEmittedClusters; };
    G4int getEnergyViolationInteraction() const { return nEnergyViolationInteraction; };

  private:
    G4int nAcceptedCollisions;
    G4int nBlockedCollisions;
    G4int nAcceptedDecays;
    G4int nBlockedDecays;
    G4double currentTime;
    G4double firstCollisionTime;
    G4double firstCollisionXSec;
    G4double firstCollisionSpectatorPosition;
    G4double firstCollisionSpectatorMomentum;
    G4bool firstCollisionIsElastic;
    std::map<AvatarType,G4int> nAvatars;
    G4int nCascading;
    G4int nEmittedClusters;
    G4int nEnergyViolationInteraction;
  };
}

#endif

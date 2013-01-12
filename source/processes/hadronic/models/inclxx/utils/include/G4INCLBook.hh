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
      nAvatars[SurfaceAvatarType] = 0;
      nAvatars[CollisionAvatarType] = 0;
      nAvatars[DecayAvatarType] = 0;
      nAvatars[ParticleEntryAvatarType] = 0;
      nCascading = 0;
    };

    void incrementAcceptedCollisions() { nAcceptedCollisions++; };
    void incrementBlockedCollisions() { nBlockedCollisions++; };
    void incrementAcceptedDecays() { nAcceptedDecays++; };
    void incrementBlockedDecays() { nBlockedDecays++; };
    void incrementAvatars(AvatarType type) { nAvatars[type]++; };
    void incrementCascading() { nCascading++; }
    void decrementCascading() { nCascading--; }

    void setFirstCollisionTime(G4double t) { firstCollisionTime = t; };
    G4double getFirstCollisionTime() { return firstCollisionTime; };

    void setFirstCollisionXSec(G4double x) { firstCollisionXSec = x; };
    G4double getFirstCollisionXSec() { return firstCollisionXSec; };

    void setCurrentTime(G4double t) { currentTime = t; };
    G4double getCurrentTime() { return currentTime; };

    G4int getAcceptedCollisions() const { return nAcceptedCollisions; };
    G4int getBlockedCollisions() const {return nBlockedCollisions; };
    G4int getAcceptedDecays() const { return nAcceptedDecays; };
    G4int getBlockedDecays() const {return nBlockedDecays; };
    G4int getAvatars(AvatarType type) const { return nAvatars.find(type)->second; };
    G4int getCascading() const { return nCascading; };

  private:
    G4int nAcceptedCollisions;
    G4int nBlockedCollisions;
    G4int nAcceptedDecays;
    G4int nBlockedDecays;
    G4double currentTime;
    G4double firstCollisionTime;
    G4double firstCollisionXSec;
    std::map<AvatarType,G4int> nAvatars;
    G4int nCascading;
  };
}

#endif

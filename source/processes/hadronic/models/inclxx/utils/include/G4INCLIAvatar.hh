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

/*
 * IAvatar.hh
 *
 *  \date 4 June 2009
 * \author Pekka Kaitaniemi
 */

#ifndef IAVATAR_HH_
#define IAVATAR_HH_

#include "G4INCLIChannel.hh"
#include "G4INCLParticle.hh"
#include "G4INCLFinalState.hh"
#include <string>

namespace G4INCL {

  enum AvatarType {SurfaceAvatarType,
		   CollisionAvatarType,
		   DecayAvatarType,
		   ParticleEntryAvatarType,
                   UnknownAvatarType};

  class IAvatar {
  public:
    IAvatar();
    IAvatar(G4double time);
    virtual ~IAvatar();

    virtual G4INCL::IChannel* getChannel() const = 0;
    G4INCL::FinalState *getFinalState();
    virtual void preInteraction() = 0;
    virtual FinalState *postInteraction(FinalState *) = 0;

    G4double getTime() const { return theTime; };

    virtual ParticleList getParticles() const = 0;

    virtual std::string dump() const = 0;

    AvatarType getType() const { return type; };
    G4bool isACollision() const { return (type==CollisionAvatarType); };
    G4bool isADecay() const { return (type==DecayAvatarType); };
    void setType(AvatarType t) { type = t; };
    long getID() const { return ID; };

    std::string toString();
  private:
    long ID;
    AvatarType type;
    static long nextID;
  protected:
    G4double theTime;
  };

  typedef std::list<IAvatar*> IAvatarList;
  typedef std::list<IAvatar*>::const_iterator IAvatarIter;
}

#endif /* IAVATAR_HH_ */

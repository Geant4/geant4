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
#include "G4INCLUnorderedVector.hh"
#include "G4INCLAllocationPool.hh"
#include <string>

#if !defined(NDEBUG) && !defined(INCLXX_IN_GEANT4_MODE)
// Force instantiation of all the std::vector<IAvatar*> methods for debugging
// purposes
namespace G4INCL {
  class IAvatar;
}
template class std::vector<G4INCL::IAvatar*>;
#endif

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

    virtual G4INCL::IChannel* getChannel() = 0;
    FinalState *getFinalState();
    void fillFinalState(FinalState *fs);
    virtual void preInteraction() = 0;
    virtual void postInteraction(FinalState *) = 0;

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
    static G4ThreadLocal long nextID;
  protected:
    G4double theTime;

    INCL_DECLARE_ALLOCATION_POOL(IAvatar)
  };

  typedef UnorderedVector<IAvatar*> IAvatarList;
  typedef UnorderedVector<IAvatar*>::const_iterator IAvatarIter;
  typedef UnorderedVector<IAvatar*>::iterator IAvatarMutableIter;

}

#endif /* IAVATAR_HH_ */

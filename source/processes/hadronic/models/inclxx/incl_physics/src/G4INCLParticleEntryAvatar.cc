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

#include "G4INCLParticleEntryAvatar.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLParticleEntryChannel.hh"

namespace G4INCL {
  //  ParticleEntryAvatar::ParticleEntryAvatar()
  //  {
  //  }

  ParticleEntryAvatar::ParticleEntryAvatar(G4double time,
					   G4INCL::Nucleus *nucleus,
					   G4INCL::Particle *particle)
    :IAvatar(time), theNucleus(nucleus), theParticle(particle)
  {
    setType(ParticleEntryAvatarType);
  }

  ParticleEntryAvatar::~ParticleEntryAvatar()
  {}

  std::string ParticleEntryAvatar::dump() const {
    std::stringstream ss;
    ss << "(avatar " << theTime <<" 'particle-entry" << std::endl
      << "(list " << std::endl
       << theParticle->dump()
      << "))" << std::endl;
    return ss.str();
  }

  FinalState* ParticleEntryAvatar::postInteraction(FinalState *fs) {
    return fs;
  }

  IChannel* ParticleEntryAvatar::getChannel() const {
    return new ParticleEntryChannel(theNucleus, theParticle);
  }
}

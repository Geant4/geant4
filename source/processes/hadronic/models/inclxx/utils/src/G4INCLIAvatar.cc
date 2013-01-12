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
 * IAvatar.cc
 *
 *  \date 4 juin 2009
 * \author Pekka Kaitaniemi
 */

#include "G4INCLIAvatar.hh"
#include <sstream>

namespace G4INCL {

  long IAvatar::nextID = 1;

  IAvatar::IAvatar() :
    type(UnknownAvatarType),
    theTime(0.)
  {
    ID = nextID;
    nextID++;
  }

  IAvatar::IAvatar(G4double time) :
    type(UnknownAvatarType),
    theTime(time)
  {
    ID = nextID;
    nextID++;
  }

  IAvatar::~IAvatar() {
  }

  std::string IAvatar::toString() {
    std::stringstream entry;
    std::stringstream particleString;
    ParticleList pl = getParticles();
    G4int numberOfParticles = 0;
    for(ParticleIter i = pl.begin(); i != pl.end(); ++i) {
      numberOfParticles++;
      particleString << (*i)->getID() << " ";
    }
    if(numberOfParticles == 1) particleString << "-1";
    entry << getID() << " "
	  << getType() << " "
	  << getTime() << " "
	  << particleString.str();
    return entry.str();
  }

  G4INCL::FinalState* IAvatar::getFinalState()
  {
    preInteraction();
    IChannel *c = getChannel();
    if( !c ) {
      return new FinalState;
    }
    FinalState *fs = c->getFinalState();
    fs = postInteraction(fs);
    delete c;
    return fs;
  }

}

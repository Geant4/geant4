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
// INCL++ revision: v5.1.2
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#ifndef G4INCLParticleEntryAvatar_hh
#define G4INCLParticleEntryAvatar_hh 1

#include "G4INCLIAvatar.hh"
#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"

namespace G4INCL {
  class ParticleEntryAvatar: public G4INCL::IAvatar {
  public:
    ParticleEntryAvatar(G4double, G4INCL::Nucleus*, G4INCL::Particle*);
    virtual ~ParticleEntryAvatar();
    virtual G4INCL::IChannel* getChannel() const;
    ParticleList getParticles() const {
      ParticleList theParticleList;
      theParticleList.push_back(theParticle);
      return theParticleList;
    };

    virtual void preInteraction() {};
    virtual FinalState *postInteraction(FinalState *);

    std::string dump() const;

    static const G4double cutNN;
    static const G4double cutNNSquared;
  private:
    Nucleus *theNucleus;
    Particle *theParticle;
  };
}

#endif

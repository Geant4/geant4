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

#include "G4INCLOmegaNElasticChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"

namespace G4INCL {

    OmegaNElasticChannel::OmegaNElasticChannel(Particle *p1, Particle *p2)
    : particle1(p1), particle2(p2)
    {

    }

    OmegaNElasticChannel::~OmegaNElasticChannel(){

    }

    void OmegaNElasticChannel::fillFinalState(FinalState *fs) {
        Particle * nucleon;
        Particle * omega;
        if(particle1->isNucleon()) {
            nucleon = particle1;
            omega = particle2;
        } else {
            nucleon = particle2;
            omega = particle1;
        }

					G4double sh=nucleon->getEnergy()+omega->getEnergy();
					G4double mn=nucleon->getMass();
					G4double me=omega->getMass();
					G4double en=(sh*sh+mn*mn-me*me)/(2*sh);
					nucleon->setEnergy(en);
					G4double ee=std::sqrt(en*en-mn*mn+me*me);
					omega->setEnergy(ee);
					G4double pn=std::sqrt(en*en-mn*mn);					

					ThreeVector mom_nucleon;

// Isotropy
					mom_nucleon = Random::normVector(pn);


					nucleon->setMomentum(mom_nucleon);
					omega->setMomentum(-mom_nucleon);

				 fs->addModifiedParticle(nucleon);
					fs->addModifiedParticle(omega);
    
				}
}

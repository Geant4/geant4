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

#include "G4INCLPiNElasticChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"

namespace G4INCL {

    PiNElasticChannel::PiNElasticChannel(Particle *p1, Particle *p2)
    : particle1(p1), particle2(p2)
    {

    }

    PiNElasticChannel::~PiNElasticChannel(){

    }

    void PiNElasticChannel::fillFinalState(FinalState *fs) {
        Particle * nucleon;
        Particle * pion;
        if(particle1->isNucleon()) {
            nucleon = particle1;
            pion = particle2;
        } else {
            nucleon = particle2;
            pion = particle1;
        }
        
        G4double bpn=8e-6;  // B-parameter of exp(Bt) - (MeV/c)^-2
        G4double px_nucleon=nucleon->getMomentum().getX();
        G4double py_nucleon=nucleon->getMomentum().getY();
        G4double pz_nucleon=nucleon->getMomentum().getZ();
        G4double cnorm1=px_nucleon*px_nucleon+py_nucleon*py_nucleon;
        G4double cnorm2=std::sqrt(cnorm1);
        G4double tnorm=cnorm1+pz_nucleon*pz_nucleon;
        G4double cnorm=std::sqrt(tnorm);
        G4double btm=std::exp(-4.*tnorm*bpn);
        G4double rndm=Random::shoot();
        G4double yrn=1.-rndm*(1.-btm);
        G4double tt=std::log(yrn)/bpn;
        G4double ctheta=1.+0.5*tt/tnorm;
        G4double stheta=std::sqrt(1.-ctheta*ctheta);
        G4double t7=1.-2.*Random::shoot();
        G4double t8=std::sqrt(1-t7*t7);
        G4double t1=-py_nucleon/cnorm2;
        G4double t2=px_nucleon/cnorm2;
        G4double t4=t2*pz_nucleon/cnorm;
        G4double t5=-t1*pz_nucleon/cnorm;
        G4double t6=-cnorm2/cnorm;
        ThreeVector mom_nucleon(
                                px_nucleon*ctheta+cnorm*stheta*(t1*t7+t4*t8),
                                py_nucleon*ctheta+cnorm*stheta*(t2*t7+t5*t8),
                                pz_nucleon*ctheta+cnorm*stheta*t6*t8
                                );
        nucleon->setMomentum(mom_nucleon);
        pion->setMomentum(-mom_nucleon);

	ParticleType startingNucleonType = nucleon->getType();
	ParticleType startingPionType = pion->getType();
	    
        G4int iso=ParticleTable::getIsospin(nucleon->getType())+ParticleTable::getIsospin(pion->getType());
        if (iso == 1 || iso == -1) {
            rndm=3*Random::shoot();
            if (rndm < 1.) {
                ParticleType nucleonType=ParticleTable::getNucleonType(-iso);
                nucleon->setType(nucleonType);
                ParticleType pionType=ParticleTable::getPionType(2*iso);
                pion->setType(pionType);
            }
            else {
                ParticleType nucleonType=ParticleTable::getNucleonType(iso);
                nucleon->setType(nucleonType);
                pion->setType(PiZero);
            }
        }
        else { // Useless
            ParticleType nucleonType=ParticleTable::getNucleonType(iso/3);
            nucleon->setType(nucleonType);
            ParticleType pionType=ParticleTable::getPionType(2*iso/3);
            pion->setType(pionType);
        }

	// Erase the parent resonance information if the nucleon or pion changes type
	if ( startingNucleonType != nucleon->getType() || startingPionType != pion->getType() ) {
	  nucleon->setParentResonancePDGCode(0);
	  nucleon->setParentResonanceID(0);
	  pion->setParentResonancePDGCode(0);
	  pion->setParentResonanceID(0);
	}
	
        fs->addModifiedParticle(nucleon);
        fs->addModifiedParticle(pion);
    }

}

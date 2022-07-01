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

#include "G4INCLPiNToEtaChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"

namespace G4INCL {

    PiNToEtaChannel::PiNToEtaChannel(Particle *p1, Particle *p2)
    : particle1(p1), particle2(p2)
    {

    }

    PiNToEtaChannel::~PiNToEtaChannel(){

    }

    void PiNToEtaChannel::fillFinalState(FinalState *fs) {
        Particle * nucleon;
        Particle * pion;
        if(particle1->isNucleon()) {
            nucleon = particle1;
            pion = particle2;
        } else {
            nucleon = particle2;
            pion = particle1;
        }

		G4int iso=ParticleTable::getIsospin(nucleon->getType())+ParticleTable::getIsospin(pion->getType());
// assert(iso == 1 || iso == -1);
		if (iso == 1) {
			nucleon->setType(Proton);
		}
		else if (iso == -1) {
			nucleon->setType(Neutron);
        }
		pion->setType(Eta);

                // Erase the parent resonance information of the nucleon and pion
	        nucleon->setParentResonancePDGCode(0);
	        nucleon->setParentResonanceID(0);
                pion->setParentResonancePDGCode(0);
                pion->setParentResonanceID(0);

		G4double sh=nucleon->getEnergy()+pion->getEnergy();
		G4double mn=nucleon->getMass();
		G4double me=pion->getMass();
		G4double en=(sh*sh+mn*mn-me*me)/(2*sh);
		nucleon->setEnergy(en);
		G4double ee=std::sqrt(en*en-mn*mn+me*me);
		pion->setEnergy(ee);
		G4double pn=std::sqrt(en*en-mn*mn);

// real distribution (from PRC 78, 025204 (2008))
 

		G4double ECM=G4INCL::KinematicsUtils::totalEnergyInCM(particle1,particle2);

		const G4double pi=std::acos(-1.0);
		G4double x1;
		G4double u1;
		G4double fteta;
		G4double teta;
		G4double fi;

		if (ECM < 1650.) {
// below 1650 MeV - angular distribution (x=cos(theta): ax^2+bx+c		
        
		G4double f1= -0.0000288627*ECM*ECM+0.09155289*ECM-72.25436;  // f(1) that is the maximum (fit on experimental data)
		G4double b1=(f1-(f1/(1.5-0.5*std::pow((ECM-1580.)/95.,2))))/2.; // ideas: 1) f(-1)=0.5f(1); 2) "power term" flattens the distribution away from ECM=1580 MeV
		G4double a1=2.5*b1; // minimum at cos(theta) = -0.2
		G4double c1=f1-3.5*b1;
				   
		G4double interg1=2.*a1/3. +2.*c1;	// (integral to normalize)   
 
		G4int passe1=0;
		while (passe1==0) {
			// Sample x from -1 to 1
			x1=Random::shoot();
			if (Random::shoot() > 0.5) x1=-x1;
			
			// Sample u from 0 to 1
			u1=Random::shoot();
			fteta=(a1*x1*x1+b1*x1+c1)/interg1;
			// The condition
			if (u1*f1/interg1 < fteta) {
				teta=std::acos(x1);
				passe1=1;
			}
		}
	}
	else {		   
// above 1650 MeV - angular distribution (x=cos(theta): (ax^2+bx+c)*(0.5+(arctan(10*(x+dev)))/pi) + vert
        
		G4double a2=-0.29;
		G4double b2=0.348;    // ax^2+bx+c: around cos(theta)=0.6 with maximum at 0.644963 (value = 0.1872666)
		G4double c2=0.0546;
		G4double dev=-0.2;  // tail close to zero from "dev" down to -1
		G4double vert=0.04; // to avoid negative differential cross sections
			
		G4double interg2=0.1716182902205207; // with the above given parameters! (integral to normalize)
		const G4double f2=1.09118088; // maximum (integral taken into account)
			
		G4int passe2=0;
		while (passe2==0) {
			// Sample x from -1 to 1
			x1=Random::shoot();
			if (Random::shoot() > 0.5) x1=-x1;
			
			// Sample u from 0 to 1
			u1=Random::shoot();
			fteta=((a2*x1*x1+b2*x1+c2)*(0.5+(std::atan(10*(x1+dev)))/pi) + vert)/interg2;
			// The condition
			if (u1*f2 < fteta) {
				teta=std::acos(x1);
				passe2=1;
			}
		}
 }
				   
		fi=(2.0*pi)*Random::shoot();		

		ThreeVector mom_nucleon(
                                pn*std::sin(teta)*std::cos(fi),
                                pn*std::sin(teta)*std::sin(fi),
                                pn*std::cos(teta)
                                );
// end real distribution			
		
		nucleon->setMomentum(-mom_nucleon);
		pion->setMomentum(mom_nucleon);
        
		fs->addModifiedParticle(nucleon);
		fs->addModifiedParticle(pion);
		}

}

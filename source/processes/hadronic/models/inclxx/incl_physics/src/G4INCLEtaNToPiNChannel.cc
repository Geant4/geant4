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

#include "G4INCLEtaNToPiNChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"

namespace G4INCL {
    
    EtaNToPiNChannel::EtaNToPiNChannel(Particle *p1, Particle *p2)
    : particle1(p1), particle2(p2)
    {
        
    }
    
    EtaNToPiNChannel::~EtaNToPiNChannel(){
        
    }
    
    void EtaNToPiNChannel::fillFinalState(FinalState *fs) {
        Particle * nucleon;
        Particle * eta;
        if(particle1->isNucleon()) {
            nucleon = particle1;
            eta = particle2;
        } else {
            nucleon = particle2;
            eta = particle1;
        }
        
        G4double plab=KinematicsUtils::momentumInLab(particle1, particle2);

        const G4double r2 = Random::shoot();
        if (nucleon->getType() == Neutron) {
            if (r2*3. < 2.) {
                nucleon->setType(Proton);
                eta->setType(PiMinus);
            }
            else {
                nucleon->setType(Neutron);
                eta->setType(PiZero);
            }
        }
        else {
            if (r2*3. < 2.) {
                nucleon->setType(Neutron);
                eta->setType(PiPlus);
            }
            else {
                nucleon->setType(Proton);
                eta->setType(PiZero);
            }
        }
        
        G4double sh=nucleon->getEnergy()+eta->getEnergy();
        G4double mn=nucleon->getMass();
        G4double me=eta->getMass();
        G4double en=(sh*sh+mn*mn-me*me)/(2*sh);
        nucleon->setEnergy(en);
        G4double ee=std::sqrt(en*en-mn*mn+me*me);
        eta->setEnergy(ee);
        G4double pn=std::sqrt(en*en-mn*mn);
                
        const G4double pi=std::acos(-1.0);
        G4double x1;
        G4double u1;
        G4double fteta;
        G4double teta;
        G4double fi;
        
        G4double a0;
        G4double a1;
        G4double a2;
        G4double a3;
        G4double a4;
        G4double a5;
        G4double a6;
        
        if (plab > 1400.) plab=1400.; // no information on angular distributions above plab=1400 MeV
        G4double p6=std::pow(plab, 6);
        G4double p5=std::pow(plab, 5);
        G4double p4=std::pow(plab, 4);
        G4double p3=std::pow(plab, 3);
        G4double p2=std::pow(plab, 2);
        G4double p1=plab;
        
        // a6
        if (plab <= 600.) {
            a6=5.721872E-18*p6 - 1.063594E-14*p5 +
            7.812226E-12*p4 - 2.947343E-09*p3 +
            5.955500E-07*p2 - 6.081534E-05*p1 + 2.418893E-03;
        }
        else {
            a6=1.549323E-18*p6 - 9.570613E-15*p5 +
            2.428560E-11*p4 - 3.237490E-08*p3 +
            2.385312E-05*p2 - 9.167580E-03*p1 + 1.426952E+00;
        }
        // a5
        if (plab <= 700.) {
            a5=-3.858406E-16*p6 + 7.397533E-13*p5 -
            5.344420E-10*p4 + 1.865842E-07*p3 -
            3.234292E-05*p2 + 2.552380E-03*p1 - 6.810842E-02;
        }
        else {
            a5=-3.775268E-17*p6 + 2.445059E-13*p5 -
            6.503137E-10*p4 + 9.065678E-07*p3 -
            6.953576E-04*p2 + 2.757524E-01*p1 - 4.328028E+01;
        }
        // a4
        if (plab <= 550.) {
            a4=-2.051840E-16*p6 + 3.858551E-13*p5 -
            3.166229E-10*p4 + 1.353545E-07*p3 -
            2.631251E-05*p2 + 2.109593E-03*p1 - 5.633076E-02;
        }
        else if (plab <= 650.) {
            a4=-1.698136E-05*p2 + 1.827203E-02*p1 - 4.482122E+00;
        }
        else {
            a4=-2.808337E-17*p6 + 1.640033E-13*p5 -
            3.820460E-10*p4 + 4.452787E-07*p3 -
            2.621981E-04*p2 + 6.530743E-02*p1 - 2.447717E+00;
        }
        // a3
        if (plab <= 700.) {
            a3=7.061866E-16*p6 - 1.356389E-12*p5 +
            9.783322E-10*p4 - 3.407333E-07*p3 +
            5.903545E-05*p2 - 4.735559E-03*p1 + 1.270435E-01;
        }
        else {
            a3=1.138088E-16*p6 - 7.459580E-13*p5 +
            2.015156E-09*p4 - 2.867416E-06*p3 +
            2.261028E-03*p2 - 9.323442E-01*p1 + 1.552846E+02;
        }
        // a2
        if (plab <= 550.) {
            a2=1.352952E-17*p6 - 3.030435E-13*p5 +
            4.624668E-10*p4 - 2.759605E-07*p3 +
            6.996373E-05*p2 - 4.745692E-03*p1 + 1.524349E-01;
        }
        else if (plab <= 700.) {
            a2=5.514651E-08*p3 - 8.734112E-05*p2 + 4.108704E-02*p1 - 5.116601E+00;
        }
        else {
            a2=5.621795E-17*p6 - 3.701960E-13*p5 +
            1.005796E-09*p4 - 1.441294E-06*p3 +
            1.146234E-03*p2 - 4.775194E-01*p1 + 8.084776E+01;
        }
        // a1
        if (plab <= 500.) {
            a1=-2.425827E-16*p6 + 4.113350E-13*p5 -
            2.342298E-10*p4 + 4.934322E-08*p3 -
            3.564530E-06*p2 + 6.516398E-04*p1 + 2.547230E-01;
        }
        else if (plab <= 700.) {
            a1=-1.824213E-10*p4 + 3.599251E-07*p3 -
            2.480862E-04*p2 + 6.894931E-02*p1 - 5.760562E+00;
        }
        else {
            a1=-5.139366E-17*p6 + 3.408224E-13*p5 -
            9.341903E-10*p4 + 1.354028E-06*p3 -
            1.093509E-03*p2 + 4.653326E-01*p1 - 8.068436E+01;
        }
        // a0
        if (plab <= 400.) {
            a0=1.160837E-13*p6 - 1.813002E-10*p5 +
            1.155391E-07*p4 - 3.862737E-05*p3 +
            7.230513E-03*p2 - 7.469799E-01*p1 + 3.830064E+01;
        }
        else if (plab <= 700.) {
            a0=2.267918E-14*p6 - 7.593899E-11*p5 +
            1.049849E-07*p4 - 7.669301E-05*p3 +
            3.123846E-02*p2 - 6.737221E+00*p1 + 6.032010E+02;
        }
        else {
            a0=-1.851188E-17*p6 + 1.281122E-13*p5 -
            3.686161E-10*p4 + 5.644116E-07*p3 -
            4.845757E-04*p2 + 2.203918E-01*p1 - 4.100383E+01;
        }
        
        G4double interg1=2.*(a6/7. + a4/5. + a2/3. + a0);	// (integral to normalize)
        G4double f1=(a6+a5+a4+a3+a2+a1+a0)/interg1; // (Max normalized)
        
        G4int passe1=0;
        while (passe1==0) {
            // Sample x from -1 to 1
            x1=Random::shoot();
            if (Random::shoot() > 0.5) x1=-x1;
            
            // Sample u from 0 to 1
            u1=Random::shoot();
            fteta=(a6*x1*x1*x1*x1*x1*x1+a5*x1*x1*x1*x1*x1+a4*x1*x1*x1*x1+a3*x1*x1*x1+a2*x1*x1+a1*x1+a0)/interg1;
            // The condition
            if (u1*f1 < fteta) {
                teta=std::acos(x1);
                //				std::cout << x1  << " " << fteta << " "<< f1/interg1 << " " << u1 << " " << interg1 << std::endl;
                passe1=1;
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
        eta->setMomentum(mom_nucleon);
        
        fs->addModifiedParticle(nucleon);
        fs->addModifiedParticle(eta);
    }
    
}

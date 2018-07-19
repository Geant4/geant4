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

#include "G4INCLEtaNElasticChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"

namespace G4INCL {

    EtaNElasticChannel::EtaNElasticChannel(Particle *p1, Particle *p2)
    : particle1(p1), particle2(p2)
    {

    }

    EtaNElasticChannel::~EtaNElasticChannel(){

    }

    void EtaNElasticChannel::fillFinalState(FinalState *fs) {
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

					G4double sh=nucleon->getEnergy()+eta->getEnergy();
					G4double mn=nucleon->getMass();
					G4double me=eta->getMass();
					G4double en=(sh*sh+mn*mn-me*me)/(2*sh);
					nucleon->setEnergy(en);
					G4double ee=std::sqrt(en*en-mn*mn+me*me);
					eta->setEnergy(ee);
					G4double pn=std::sqrt(en*en-mn*mn);					

					ThreeVector mom_nucleon;

					if (plab < 250.) {
// Isotropy
					mom_nucleon = Random::normVector(pn);
					}					

// From Kamano
     else {

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
						if (plab < 300.) {
       a6=-8.384000E-08*p1 - 1.15452E-04;
						}
      else if (plab < 500.){
							a6=1.593966E-13*p4 - 2.619560E-10*p3 + 1.564701E-07*p2 - 3.986627E-05*p1 + 3.622575E-03;
      }
						else {
							a6=6.143615E-20*p6 - 3.157181E-16*p5 + 6.348289E-13*p4 - 6.117961E-10*p3 + 2.764542E-07*p2 - 4.391048E-05*p1 - 1.443857E-03;
						}
// a5
						if (plab < 650.) {
							a5=-9.021076E-18*p6 + 2.176771E-14*p5 - 2.136095E-11*p4 + 1.100580E-08*p3 - 3.150857E-06*p2 + 4.761016E-04*p1 - 2.969608E-02;
						}
      else if (plab < 950.){
							a5=4.424756E-18*p6 - 1.756295E-14*p5 + 2.625428E-11*p4 - 1.678272E-08*p3 + 2.227237E-06*p2 + 2.146666E-03*p1 - 7.065712E-01;
      }
						else {
							a5=2.209585E-19*p6 - 1.546647E-15*p5 + 4.578142E-12*p4 - 7.303856E-09*p3 + 6.604074E-06*p2 - 3.205628E-03*p1 + 6.534893E-01;
						}
// a4
						if (plab < 700.) {
							a4=4.826684E-17*p6 - 1.534471E-13*p5 + 1.907868E-10*p4 - 1.192317E-07*p3 + 3.988902E-05*p2 - 6.822100E-03*p1 + 4.684685E-01;
						}
						else {
							a4=-3.245143E-18*p6 + 2.174395E-14*p5 - 6.012288E-11*p4 + 8.772790E-08*p3 - 7.113554E-05*p2 + 3.029285E-02*p1 - 5.237677E+00;
						}
// a3
						if (plab < 650.) {
							a3=3.783071E-17*p6 - 1.151454E-13*p5 + 1.357165E-10*p4 - 8.036891E-08*p3 + 2.572396E-05*p2 - 4.245566E-03*p1 + 2.832772E-01;
						}
						else {
							a3=-5.063316E-18*p6 + 3.223757E-14*p5 - 8.435635E-11*p4 + 1.159487E-07*p3 - 8.812510E-05*p2 + 3.500692E-02*p1 - 5.624556E+00;
						}
// a2
						if (plab < 500.) {
							a2=-6.085067E-14*p5 + 1.354078E-10*p4 - 1.124158E-07*p3 + 4.292106E-05*p2 - 7.218145E-03*p1 + 4.584962E-01;
						}
						else if (plab < 750.) {
							a2= 9.512730E-11*p4 - 2.362724E-07*p3 + 2.171883E-04*p2 - 8.742722E-02*p1 + 1.309433E+01;
						}
						else {
							a2=-4.228889E-18*p6 + 2.798222E-14*p5 - 7.640831E-11*p4 + 1.100124E-07*p3 - 8.778573E-05*p2 + 3.652772E-02*p1 - 6.025497E+00;
						}
// a1
						if (plab < 500.) {
							a1=-1.524408E-14*p5 + 3.007021E-11*p4 - 2.129570E-08*p3 + 5.607250E-06*p2 - 3.001598E-04*p1 + 8.701280E-04;
						}
						else if (plab < 750.) {
							a1=-3.255396E-11*p4 + 8.168681E-08*p3 - 7.447474E-05*p2 + 2.917630E-02*p1 - 4.152037E+00;
						}
						else {
							a1=9.964504E-19*p6 - 6.380168E-15*p5 + 1.638691E-11*p4 - 2.107063E-08*p3 + 1.347462E-05*p2 - 3.318304E-03*p1 - 5.030932E-02;
						}
// a0
							a0=-3.220143E-17*p6 + 1.789654E-13*p5 - 3.912863E-10*p4 + 4.181510E-07*p3 - 2.147259E-04*p2 + 3.856266E-02*p1 + 2.609971E+00;
						
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
						
						ThreeVector mom_nucleon1(
																														pn*std::sin(teta)*std::cos(fi),
																														pn*std::sin(teta)*std::sin(fi),
																														pn*std::cos(teta)
																														);
						
						mom_nucleon = -mom_nucleon1	;
					
					}

					nucleon->setMomentum(mom_nucleon);
					eta->setMomentum(-mom_nucleon);

				 fs->addModifiedParticle(nucleon);
					fs->addModifiedParticle(eta);
    
				}
}

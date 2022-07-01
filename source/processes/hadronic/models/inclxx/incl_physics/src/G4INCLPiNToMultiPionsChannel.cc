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

#include "G4INCLPiNToMultiPionsChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {

  const G4double PiNToMultiPionsChannel::angularSlope = 15.;

  PiNToMultiPionsChannel::PiNToMultiPionsChannel(const G4int npi, Particle *p1, Particle *p2)
    : npion(npi),
    ind2(0),
    particle1(p1),
    particle2(p2)
  {
    std::fill(isosp, isosp+4, 0);
  }

  PiNToMultiPionsChannel::~PiNToMultiPionsChannel(){

  }

  void PiNToMultiPionsChannel::fillFinalState(FinalState *fs) {

// assert(npion > 1 && npion < 5);

    Particle * nucleon;
    Particle * pion;
    if(particle1->isNucleon()) {
      nucleon = particle1;
      pion = particle2;
    } else {
      nucleon = particle2;
      pion = particle1;
    }

      // Erase the parent resonance information of the nucleon and pion
      nucleon->setParentResonancePDGCode(0);
      nucleon->setParentResonanceID(0);
      pion->setParentResonancePDGCode(0);
      pion->setParentResonanceID(0);
    
      G4int ipi=ParticleTable::getIsospin(pion->getType());
      ind2=ParticleTable::getIsospin(nucleon->getType());

      ParticleList list;
      list.push_back(nucleon);
      list.push_back(pion);
      fs->addModifiedParticle(nucleon);
      fs->addModifiedParticle(pion);

      isospinRepartition(ipi);

      const ParticleType tn=ParticleTable::getNucleonType(ind2);
      nucleon->setType(tn);
      ParticleType pionType=ParticleTable::getPionType(isosp[0]);
      pion->setType(pionType);
      const ThreeVector &rcolpion = pion->getPosition();
      const ThreeVector zero;
      for(G4int i=1; i<npion; ++i) {
        pionType=ParticleTable::getPionType(isosp[i]);
        Particle *newPion = new Particle(pionType,zero,rcolpion);
        newPion->setType(pionType);
        list.push_back(newPion);
        fs->addCreatedParticle(newPion);
      }

      const G4double sqrtS = KinematicsUtils::totalEnergyInCM(nucleon, pion);
      PhaseSpaceGenerator::generateBiased(sqrtS, list, 0, angularSlope);

  }

    void PiNToMultiPionsChannel::isospinRepartition(G4int ipi) {
        G4double rjcd=Random::shoot();
        const G4int itot=ipi*ind2;

        isosp[1]=ipi;
        if (npion != 3) {
            if (npion == 4){
              const G4double r2 = Random::shoot();
                if (r2*3. > 2.) {
                    isosp[2]= 0;
                    isosp[3]= 0;
                }
                else {
                    isosp[2]=  2;
                    isosp[3]= -2;
                }
            }

            if (itot == 2) {
                // CAS PI+ P ET PI- n
                rjcd *= 5.;
                if (rjcd > 3.) {
                    //      PI+ PI+ N ET PI- PI- P
                    isosp[0]=2*ind2;
                    isosp[1]=ipi;
                    ind2=-ind2;
                }
                else {
                    //      PI+ PI0 P ET PI- PI0 N
                    isosp[0]=0;
                    isosp[1]=ipi;
                }
            }
            else if (itot == 0) {
                // CAS PI0 P ET PI0 N
                rjcd *= 90.;
                if (rjcd > 13.) {
                    if (rjcd > 52.) {
                        //     PI+ PI0 N ET PI- PI0 P
                        isosp[0]=2*ind2;
                        isosp[1]=0;
                        ind2=-ind2;
                    }
                    else {
                        //     PI+ PI- P ET PI+ PI- N
                        isosp[1]=-2;
                        isosp[0]=2;
                    }
                }
                else {
                    //     PI0 PI0 P ET PI0 PI0 N
                    isosp[0]=0;
                    isosp[1]=0;
                }
            }
            else if (itot == -2) {
                // CAS PI- P ET PI+ N
                rjcd *= 45.;
                if (rjcd > 17.) {
                    if (rjcd > 24.) {
                        //     PI+ PI- N ET PI+ PI- P
                        isosp[0]=2*ind2;
                        ind2=-ind2;
                    }
                    else {
                        //     PI0 PI0 N ET PI0 PI0 P
                        isosp[0]=0;
                        isosp[1]=0;
                        ind2=-ind2;
                    }
                }
                else
                    //     PI- PI0 P ET PI+ PI0 N
                    isosp[0]=0;
            }
        } // if (npion != 3)
        else {
            // PI N -> PI PI PI N
            //          IF (IPI*IND2) 20,21,22
            if (itot == -2) {
                // CAS PI- P ET PI+ N
                rjcd *= 135.;
                if (rjcd <= 28.) {
                    //     PI0 PI0 PI0 N ET PI0 PI0 PI0 P
                    isosp[0]=0;
                    isosp[1]=0;
                    isosp[2]=0;
                    ind2=-ind2;
                }
                else {
                    if (rjcd <= 84.) {
                        //     PI+ PI- PI0 N ET PI- PI+ PI0 P
                        isosp[0]=2*ind2;
                        isosp[2]=0;
                        ind2=-ind2;
                    }
                    else {
                        if (rjcd <= 118.) {
                            //     PI- PI- PI+ P ET PI- PI+ PI+ N
                            isosp[0]=ipi;
                            isosp[2]=-ipi;
                        }
                        else {
                            //     PI- PI0 PI0 P ET PI0 PI0 PI+ N
                            isosp[0]=0;
                            isosp[2]=0;
                        }
                    }
                }
            }
            else if (itot == 0) {
                // CAS PI0 P ET PI0 N
                rjcd *= 270.;
                if (rjcd <= 39.) {
                    //     PI0 PI0 PI0 P ET PI0 PI0 PI0 N
                    isosp[0]=0;
                    isosp[2]=0;
                }
                else {
                    if (rjcd <= 156.) {
                        //     PI+ PI- PI0 P ET PI- PI+ PI0 N
                        isosp[0]=2;
                        isosp[2]=-2;
                    }
                    else {
                        if (rjcd <= 194.) {
                            //     PI+ PI0 PI0 N ET PI0 PI0 PI- P
                            isosp[0]=0;
                            isosp[2]=2*ind2;
                            ind2=-ind2;
                        }
                        else {
                            //     PI- PI+ PI+ N ET PI- PI- PI+ P
                            isosp[0]=2*ind2;
                            isosp[1]=2*ind2;
                            isosp[2]=-2*ind2;
                            ind2=-ind2;
                        }
                    }
                }
            }
            else if (itot == 2) {
                // CAS PI+ P ET PI- N
                rjcd *= 5.;
                if (rjcd <= 2.) {
                    //     PI+ P PI0 PI0 ET PI- N PI0 PI0
                    isosp[0]=0;
                    isosp[2]=0;
                }
                else {
                    if (rjcd <= 3.) {
                        //     PI+ P PI+ PI- ET PI- N PI+ PI-
                        isosp[0]=-2;
                        isosp[2]=2;
                    }
                    else {
                        //     PI+ N PI+ PI0 ET PI- P PI0 PI-
                        isosp[0]=2*ind2;
                        isosp[2]=0;
                        ind2=-ind2;
                    }
                }
            }
        }

        std::shuffle(isosp,isosp+npion,Random::getAdapter()); // isospin randomly distributed
    }

}

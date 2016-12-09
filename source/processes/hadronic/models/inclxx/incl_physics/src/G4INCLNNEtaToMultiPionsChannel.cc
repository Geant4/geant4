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

#include "G4INCLNNEtaToMultiPionsChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {

  const G4double NNEtaToMultiPionsChannel::angularSlope = 6.;

  NNEtaToMultiPionsChannel::NNEtaToMultiPionsChannel(const G4int npi, Particle *p1, Particle *p2)
    : npion(npi),
    iso1(0),
    iso2(0),
    particle1(p1),
    particle2(p2)
  {
    std::fill(isosp, isosp+4, 0);
  }

  NNEtaToMultiPionsChannel::~NNEtaToMultiPionsChannel(){

  }

  void NNEtaToMultiPionsChannel::fillFinalState(FinalState *fs) {
// assert(npion > 0 && npion < 5);

      iso1=ParticleTable::getIsospin(particle1->getType());
      iso2=ParticleTable::getIsospin(particle2->getType());

      ParticleList list;
      list.push_back(particle1);
      list.push_back(particle2);
      fs->addModifiedParticle(particle1);
      fs->addModifiedParticle(particle2);

      isospinRepartition();

      const ParticleType tn1=ParticleTable::getNucleonType(iso1);
      particle1->setType(tn1);
      const ParticleType tn2=ParticleTable::getNucleonType(iso2);
      particle2->setType(tn2);
      const ThreeVector &rcolnucleon1 = particle1->getPosition();
      const ThreeVector &rcolnucleon2 = particle2->getPosition();
      const ThreeVector rcol = (rcolnucleon1+rcolnucleon2)*0.5;
      const ThreeVector zero;
      for(G4int i=0; i<npion; ++i) {
        const ParticleType pionType=ParticleTable::getPionType(isosp[i]);
        Particle *pion = new Particle(pionType,zero,rcol);
        list.push_back(pion);
        fs->addCreatedParticle(pion);
      }
      Particle *eta = new Particle(Eta,zero,rcol);
      list.push_back(eta);
      fs->addCreatedParticle(eta);
   
      const G4double sqrtS = KinematicsUtils::totalEnergyInCM(particle1, particle2);
      G4int biasIndex = ((Random::shoot()<0.5) ? 0 : 1);
      PhaseSpaceGenerator::generateBiased(sqrtS, list, biasIndex, angularSlope);

  }

    void NNEtaToMultiPionsChannel::isospinRepartition() {
        const G4double rjcd=Random::shoot();
        G4double p;
        const G4int itot=iso1+iso2;

        if (npion == 1) {
            p=3.*rjcd;
            if (p < 1.)      pn_ppPim();
            else if (p < 2.) pn_pnPi0();
            else            pn_nnPip();
        }
        else if (npion == 2) {
            if (itot == 2) {
                p=20.*rjcd;
                if (p >= 14.)      pp_nnPipPip();
                else if (p >= 11.) pp_pnPipPi0();
                else if (p >= 7.)  pp_ppPi0Pi0();
                else              pp_ppPipPim();
            }
            else if (itot == -2) {
                p=20.*rjcd;
                if (p >= 14.)      nn_ppPimPim();
                else if (p >= 11.) nn_pnPimPi0();
                else if (p >= 7.)  nn_nnPi0Pi0();
                else              nn_nnPipPim();
            }
            else  {
                G4double pp=Random::shoot();
                if (pp > 0.5) {
                    p=3.*rjcd;
                    if (p < 2.) {
                        pn_pnPipPim();
                    }
                    else {
                        pn_pnPi0Pi0();
                    }
                }
                else {
                    p=60.*rjcd;
                    if (p >= 51.)      pn_nnPipPi0();
                    else if (p >= 33.) pn_pnPi0Pi0();
                    else if (p >= 9.)  pn_pnPipPim();
                    else              pn_ppPimPi0();
                }
            }
        }
        else if (npion == 3) {
            p=60.*rjcd;
            if (itot == 2) {
                if (p >= 42.)      pp_nnPipPipPi0();
                else if (p >= 39.) pp_pnPipPi0Pi0();
                else if (p >= 33.) pp_pnPipPipPim();
                else if (p >= 22.) pp_ppPi0Pi0Pi0();
                else              pp_ppPipPimPi0();
            }
            else if (itot == -2) {
                if (p >= 42.)      nn_ppPimPimPi0();
                else if (p >= 39.) nn_pnPimPi0Pi0();
                else if (p >= 33.) nn_pnPipPimPim();
                else if (p >= 22.) nn_nnPi0Pi0Pi0();
                else              nn_nnPipPimPi0();
            }
            else {
                if (p >= 57.)      pn_nnPipPi0Pi0();
                else if (p >= 51.) pn_nnPipPipPim();
                else if (p >= 37.) pn_pnPi0Pi0Pi0();
                else if (p >= 9.)  pn_pnPi0PipPim();
                else if (p >= 6.)  pn_ppPimPi0Pi0();
                else              pn_ppPimPimPip();

            }
        }
        else if (npion == 4) {
            p=60.*rjcd;
            if (itot == 2) {
                if (p >= 48.)      pp_nnPipPipPipPim();
                else if (p >= 42.) pp_nnPipPipPi0Pi0();
                else if (p >= 36.) pp_pnPipPipPi0Pim();
                else if (p >= 33.) pp_pnPipPi0Pi0Pi0();
                else if (p >= 19.) pp_ppPipPipPimPim();
                else if (p >= 4.)  pp_ppPipPi0Pi0Pim();
                else              pp_ppPi0Pi0Pi0Pi0();
            }
            else if (itot == -2) {
                if (p >= 48.)      nn_ppPipPimPimPim();
                else if (p >= 42.) nn_ppPi0Pi0PimPim();
                else if (p >= 36.) nn_pnPipPi0PimPim();
                else if (p >= 33.) nn_pnPi0Pi0Pi0Pim();
                else if (p >= 19.) nn_nnPipPipPimPim();
                else if (p >= 4.)  nn_nnPipPi0Pi0Pim();
                else              nn_nnPi0Pi0Pi0Pi0();
            }
            else {
                G4double pp=Random::shoot();
                if (pp > 0.5) {
                    p=9.*rjcd;
                    if (p < 1.)      pn_pnPi0Pi0Pi0Pi0();
                    else if (p < 5.) pn_pnPipPi0Pi0Pim();
                    else            pn_pnPipPipPimPim();
                }
                else {
                    if (p < 3.)       pn_ppPi0Pi0Pi0Pim();
                    else if (p < 9.)  pn_ppPipPi0PimPim();
                    else if (p < 15.) pn_pnPi0Pi0Pi0Pi0();
                    else if (p < 35.) pn_pnPipPi0Pi0Pim();
                    else if (p < 51.) pn_pnPipPipPimPim();
                    else if (p < 54.) pn_nnPipPi0Pi0Pi0();
                    else             pn_nnPipPipPi0Pim();
                }
            }
        }

        std::random_shuffle(isosp,isosp+npion,Random::getAdapter());
        inter2Part(0.5);
    }


    void NNEtaToMultiPionsChannel::pn_ppPim() {
        isosp[0]=-2;
        iso1=1;
        iso2=1;
    }
    void NNEtaToMultiPionsChannel::pn_pnPi0() {
        isosp[0]=0;
    }
    void NNEtaToMultiPionsChannel::pn_nnPip() {
        isosp[0]=2;
        iso1=-1;
        iso2=-1;
    }
    void NNEtaToMultiPionsChannel::pp_nnPipPip() {
        isosp[0]=2;
        isosp[1]=2;
        iso1=-1;
        iso2=-1;
    }
    void NNEtaToMultiPionsChannel::nn_ppPimPim() {
        isosp[0]=-2;
        isosp[1]=-2;
        iso1=1;
        iso2=1;
    }
    void NNEtaToMultiPionsChannel::pn_pnPipPim() {
        isosp[0]=2;
        isosp[1]=-2;
    }
    void NNEtaToMultiPionsChannel::pn_pnPi0Pi0() {
        isosp[0]=0;
        isosp[1]=0;
    }
    void NNEtaToMultiPionsChannel::pp_ppPipPim() {
        isosp[0]=2;
        isosp[1]=-2;
    }
    void NNEtaToMultiPionsChannel::nn_nnPipPim() {
        isosp[0]=2;
        isosp[1]=-2;
    }
    void NNEtaToMultiPionsChannel::pp_ppPi0Pi0() {
        isosp[0]=0;
        isosp[1]=0;
    }
    void NNEtaToMultiPionsChannel::nn_nnPi0Pi0() {
        isosp[0]=0;
        isosp[1]=0;
    }
    void NNEtaToMultiPionsChannel::pp_pnPipPi0() {
        isosp[0]=2;
        isosp[1]=0;
        iso1=1;
        iso2=-1;
    }
    void NNEtaToMultiPionsChannel::pn_ppPimPi0() {
        isosp[0]=-2;
        isosp[1]=0;
        iso1=1;
        iso2=1;
    }
    void NNEtaToMultiPionsChannel::pn_nnPipPi0() {
        isosp[0]=2;
        isosp[1]=0;
        iso1=-1;
        iso2=-1;
    }
    void NNEtaToMultiPionsChannel::nn_pnPimPi0() {
        isosp[0]=-2;
        isosp[1]=0;
        iso1=1;
        iso2=-1;
    }
    void NNEtaToMultiPionsChannel::pp_pnPipPi0Pi0() {
        isosp[0]=2;
        isosp[1]=0;
        isosp[2]=0;
        iso1=1;
        iso2=-1;
    }
    void NNEtaToMultiPionsChannel::nn_pnPimPi0Pi0() {
        isosp[0]=-2;
        isosp[1]=0;
        isosp[2]=0;
        iso1=1;
        iso2=-1;
    }
    void NNEtaToMultiPionsChannel::pn_nnPipPi0Pi0() {
        isosp[0]=2;
        isosp[1]=0;
        isosp[2]=0;
        iso1=-1;
        iso2=-1;
    }
    void NNEtaToMultiPionsChannel::pp_ppPipPimPi0() {
        isosp[0]=2;
        isosp[1]=-2;
        isosp[2]=0;
    }
    void NNEtaToMultiPionsChannel::nn_nnPipPimPi0() {
        isosp[0]=2;
        isosp[1]=-2;
        isosp[2]=0;
    }
    void NNEtaToMultiPionsChannel::pp_ppPi0Pi0Pi0() {
        isosp[0]=0;
        isosp[1]=0;
        isosp[2]=0;
    }
    void NNEtaToMultiPionsChannel::nn_nnPi0Pi0Pi0() {
        isosp[0]=0;
        isosp[1]=0;
        isosp[2]=0;
    }
    void NNEtaToMultiPionsChannel::pp_pnPipPipPim() {
        isosp[0]=2;
        isosp[1]=2;
        isosp[2]=-2;
        iso1=1;
        iso2=-1;
    }
    void NNEtaToMultiPionsChannel::pp_nnPipPipPi0() {
        isosp[0]=2;
        isosp[1]=2;
        isosp[2]=0;
        iso1=-1;
        iso2=-1;
    }
    void NNEtaToMultiPionsChannel::pn_ppPimPi0Pi0() {
        isosp[0]=-2;
        isosp[1]=0;
        isosp[2]=0;
        iso1=1;
        iso2=1;
    }
    void NNEtaToMultiPionsChannel::pn_ppPimPimPip() {
        isosp[0]=-2;
        isosp[1]=-2;
        isosp[2]=2;
        iso1=1;
        iso2=1;
    }
    void NNEtaToMultiPionsChannel::pn_pnPi0PipPim() {
        isosp[0]=0;
        isosp[1]=2;
        isosp[2]=-2;
    }
    void NNEtaToMultiPionsChannel::pn_pnPi0Pi0Pi0() {
        isosp[0]=0;
        isosp[1]=0;
        isosp[2]=0;
    }
    void NNEtaToMultiPionsChannel::pn_nnPipPipPim() {
        isosp[0]=2;
        isosp[1]=2;
        isosp[2]=-2;
        iso1=-1;
        iso2=-1;
    }
    void NNEtaToMultiPionsChannel::nn_pnPipPimPim() {
        isosp[0]=2;
        isosp[1]=-2;
        isosp[2]=-2;
        iso1=1;
        iso2=-1;
    }
    void NNEtaToMultiPionsChannel::nn_ppPimPimPi0() {
        isosp[0]=-2;
        isosp[1]=-2;
        isosp[2]=0;
        iso1=1;
        iso2=1;
    }
    void NNEtaToMultiPionsChannel::pp_nnPipPipPi0Pi0() {
        isosp[0]=2;
        isosp[1]=2;
        isosp[2]=0;
        isosp[3]=0;
        iso1=-1;
        iso2=-1;
    }
    void NNEtaToMultiPionsChannel::pp_nnPipPipPipPim() {
        isosp[0]=2;
        isosp[1]=2;
        isosp[2]=2;
        isosp[3]=-2;
        iso1=-1;
        iso2=-1;
    }
    void NNEtaToMultiPionsChannel::nn_ppPi0Pi0PimPim() {
        isosp[0]=0;
        isosp[1]=0;
        isosp[2]=-2;
        isosp[3]=-2;
        iso1=1;
        iso2=1;
    }
    void NNEtaToMultiPionsChannel::nn_ppPipPimPimPim() {
        isosp[0]=2;
        isosp[1]=-2;
        isosp[2]=-2;
        isosp[3]=-2;
        iso1=1;
        iso2=1;
    }
    void NNEtaToMultiPionsChannel::pp_ppPi0Pi0Pi0Pi0() {
        isosp[0]=0;
        isosp[1]=0;
        isosp[2]=0;
        isosp[3]=0;
    }
    void NNEtaToMultiPionsChannel::nn_nnPi0Pi0Pi0Pi0() {
        isosp[0]=0;
        isosp[1]=0;
        isosp[2]=0;
        isosp[3]=0;
    }
    void NNEtaToMultiPionsChannel::pn_pnPi0Pi0Pi0Pi0() {
        isosp[0]=0;
        isosp[1]=0;
        isosp[2]=0;
        isosp[3]=0;
    }
    void NNEtaToMultiPionsChannel::pp_ppPipPi0Pi0Pim() {
        isosp[0]=2;
        isosp[1]=0;
        isosp[2]=0;
        isosp[3]=-2;
    }
    void NNEtaToMultiPionsChannel::nn_nnPipPi0Pi0Pim() {
        isosp[0]=2;
        isosp[1]=0;
        isosp[2]=0;
        isosp[3]=-2;
    }
    void NNEtaToMultiPionsChannel::pn_pnPipPi0Pi0Pim() {
        isosp[0]=2;
        isosp[1]=0;
        isosp[2]=0;
        isosp[3]=-2;
    }
    void NNEtaToMultiPionsChannel::pp_ppPipPipPimPim() {
        isosp[0]=2;
        isosp[1]=2;
        isosp[2]=-2;
        isosp[3]=-2;
    }
    void NNEtaToMultiPionsChannel::nn_nnPipPipPimPim() {
        isosp[0]=2;
        isosp[1]=2;
        isosp[2]=-2;
        isosp[3]=-2;
    }
    void NNEtaToMultiPionsChannel::pn_pnPipPipPimPim() {
        isosp[0]=2;
        isosp[1]=2;
        isosp[2]=-2;
        isosp[3]=-2;
    }
    void NNEtaToMultiPionsChannel::pp_pnPipPi0Pi0Pi0() {
        isosp[0]=2;
        isosp[1]=0;
        isosp[2]=0;
        isosp[3]=0;
        iso1=1;
        iso2=-1;
    }
    void NNEtaToMultiPionsChannel::pn_nnPipPi0Pi0Pi0() {
        isosp[0]=2;
        isosp[1]=0;
        isosp[2]=0;
        isosp[3]=0;
        iso1=-1;
        iso2=-1;
    }
    void NNEtaToMultiPionsChannel::pp_nnPipPi0Pi0Pi0() {
        isosp[0]=2;
        isosp[1]=0;
        isosp[2]=0;
        isosp[3]=0;
        iso1=-1;
        iso2=-1;
    }
    void NNEtaToMultiPionsChannel::pp_pnPipPipPi0Pim() {
        isosp[0]=2;
        isosp[1]=2;
        isosp[2]=0;
        isosp[3]=-2;
        iso1=1;
        iso2=-1;
    }
    void NNEtaToMultiPionsChannel::pn_nnPipPipPi0Pim() {
        isosp[0]=2;
        isosp[1]=2;
        isosp[2]=0;
        isosp[3]=-2;
        iso1=-1;
        iso2=-1;
    }
    void NNEtaToMultiPionsChannel::pp_nnPipPipPi0Pim() {
        isosp[0]=2;
        isosp[1]=2;
        isosp[2]=0;
        isosp[3]=-2;
        iso1=-1;
        iso2=-1;
    }
    void NNEtaToMultiPionsChannel::nn_pnPi0Pi0Pi0Pim() {
        isosp[0]=0;
        isosp[1]=0;
        isosp[2]=0;
        isosp[3]=-2;
        iso1=1;
        iso2=-1;
    }
    void NNEtaToMultiPionsChannel::pn_ppPi0Pi0Pi0Pim() {
        isosp[0]=0;
        isosp[1]=0;
        isosp[2]=0;
        isosp[3]=-2;
        iso1=1;
        iso2=1;
    }
    void NNEtaToMultiPionsChannel::nn_pnPipPi0PimPim() {
        isosp[0]=2;
        isosp[1]=0;
        isosp[2]=-2;
        isosp[3]=-2;
        iso1=1;
        iso2=-1;
    }
    void NNEtaToMultiPionsChannel::pn_ppPipPi0PimPim() {
        isosp[0]=2;
        isosp[1]=0;
        isosp[2]=-2;
        isosp[3]=-2;
        iso1=1;
        iso2=1;
    }

    void NNEtaToMultiPionsChannel::inter2Part(const G4double p) {

        if (Random::shoot() < p) std::swap(iso1,iso2);

    }


}

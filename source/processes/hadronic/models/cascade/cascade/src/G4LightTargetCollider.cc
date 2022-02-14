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
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   G4LightTargetCollider.cc                                          //
//  Date:   30 September 2019                                                 //
//  Author: Dennis Wright (SLAC)                                              //
//                                                                            //
//  Description: model for collision of elementary particles with light       //
//               targets (H, D, T, 3He)                                       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4LightTargetCollider.hh"
#include "G4CascadeChannel.hh"
#include "G4CascadeChannelTables.hh"
#include "G4CascadeCheckBalance.hh"
#include "G4CollisionOutput.hh"
#include "G4ElementaryParticleCollider.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4NucleiModel.hh"
#include "G4LorentzConvertor.hh"
#include "G4Deuteron.hh"
#include "G4Gamma.hh"
#include "G4PionZero.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4RandomDirection.hh"


G4LightTargetCollider::G4LightTargetCollider()
 : G4CascadeColliderBase("G4LightTargetCollider"),
   theElementaryParticleCollider(new G4ElementaryParticleCollider)
{
  mP = G4Proton::Proton()->GetPDGMass()/CLHEP::GeV;
  mN = G4Neutron::Neutron()->GetPDGMass()/CLHEP::GeV;
  mD = G4Deuteron::Deuteron()->GetPDGMass()/CLHEP::GeV;
  pFermiD = 0.045;  // Fermi momentum of nucleon in deuteron Hulthen potential
}

G4LightTargetCollider::~G4LightTargetCollider() {
  delete theElementaryParticleCollider;
}


// Set verbosity and pass on to member objects
void G4LightTargetCollider::setVerboseLevel(G4int verbose) {
  G4CascadeColliderBase::setVerboseLevel(verbose);
  theElementaryParticleCollider->setVerboseLevel(verboseLevel);
  output.setVerboseLevel(verboseLevel);
}


void G4LightTargetCollider::collide(G4InuclParticle* bullet,
                                    G4InuclParticle* target,
                                    G4CollisionOutput& globalOutput)
{
  if (verboseLevel) {
    G4cout << " >>> G4LightTargetCollider::collide" << G4endl;
    G4cout << "     Projectile: " << bullet->getDefinition()->GetParticleName() << G4endl; 
    G4cout << "     Target: " << target->getDefinition()->GetParticleName() << G4endl;
  }

  G4double ke = bullet->getKineticEnergy();

  if (target->getDefinition() == G4Proton::Proton() ) {
    if (ke < 0.1447) {
      // Below threshold lab energy for pi0 creation 
      globalOutput.trivialise(bullet, target);
    } else {
      // Need inelastic cross section in this class if we want to replace ElementaryParticleCollider 
      // with SingleNucleonScattering
      theElementaryParticleCollider->collide(bullet, target, globalOutput);
      if (globalOutput.numberOfOutgoingParticles() == 0) globalOutput.trivialise(bullet, target);
    }

  } else if (target->getDefinition() == G4Deuteron::Deuteron() ) {

    if (ke < mP + mN - mD) {
      // Should not happen as long as inelastic cross section is zero
      G4Exception("G4LightTargetCollider::collide()","HAD_BERT_201",
                  JustWarning, "Projectile energy below reaction threshold");
      globalOutput.trivialise(bullet, target);

    } else {
      // Get p, n and deuteron cross sections; use lab energy to access 
      G4double gammaPXS = G4CascadeChannelTables::GetTable(9)->getCrossSection(ke);
      G4double gammaNXS = G4CascadeChannelTables::GetTable(18)->getCrossSection(ke);
      G4double gammaDXS = GammaDCrossSection(ke);

      G4double probP = 0.0;
      G4double probN = 0.0;
      // Highest threshold is 0.152 (for gamma p -> n pi+)
      // Because of Fermi momentum in deuteron, raise this to 0.159
      if (ke > 0.159) {
        G4double totalDXS = gammaPXS + gammaNXS + gammaDXS;
        probP = gammaPXS/totalDXS;
        probN = (gammaPXS+gammaNXS)/totalDXS;
      }

      G4double rndm = G4UniformRand();
      if (rndm < probP) {
        // Generate Fermi momenta of bullet and target
        G4ThreeVector fermiMomentum = pFermiD*G4RandomDirection();
        G4LorentzVector protonMomentum(fermiMomentum, std::sqrt(mP*mP + pFermiD*pFermiD) );
        G4LorentzVector neutronMomentum(-fermiMomentum, std::sqrt(mN*mN + pFermiD*pFermiD) );

        G4LorentzVector bulletMomentum = bullet->getMomentum();
        G4ThreeVector betacm = bulletMomentum.findBoostToCM(protonMomentum);

        // First boost bullet and target so that target is at rest
        G4ThreeVector toProtonRest = -protonMomentum.boostVector();
        protonMomentum.boost(toProtonRest);
        bulletMomentum.boost(toProtonRest);

        G4InuclElementaryParticle projectile(bulletMomentum, bullet->getDefinition() );
        G4InuclElementaryParticle targetNucleon(protonMomentum, G4Proton::Proton() );
        G4InuclElementaryParticle spectatorNucleon(neutronMomentum, G4Neutron::Neutron() );
        ScatteringProducts products = SingleNucleonScattering(projectile, targetNucleon);

        // Particles from SingleNucleonScattering are in CM frame of projectile 
        // and moving proton.  Transform back to lab frame with -betacm, then
        // add them to outgoing list.
        globalOutput.reset();
        G4LorentzVector temp;
        for (G4int i = 0; i < G4int(products.size()); i++) {
          temp = products[i].getMomentum();
          temp.boost(-betacm);
          products[i].setMomentum(temp);
          globalOutput.addOutgoingParticle(products[i]);
        }

        // Add the recoil nucleon unmodified
        globalOutput.addOutgoingParticle(spectatorNucleon);

      } else if (rndm < probN) {
        G4ThreeVector fermiMomentum = pFermiD*G4RandomDirection();
        G4LorentzVector protonMomentum(fermiMomentum, std::sqrt(mP*mP + pFermiD*pFermiD) );
        G4LorentzVector neutronMomentum(-fermiMomentum, std::sqrt(mN*mN + pFermiD*pFermiD) );

        G4LorentzVector bulletMomentum = bullet->getMomentum();
        G4ThreeVector betacm = bulletMomentum.findBoostToCM(neutronMomentum);

        // First boost bullet and target so that target is at rest
        G4ThreeVector toNeutronRest = -neutronMomentum.boostVector();
        neutronMomentum.boost(toNeutronRest);
        bulletMomentum.boost(toNeutronRest);

        G4InuclElementaryParticle projectile(bulletMomentum, bullet->getDefinition() );
        G4InuclElementaryParticle targetNucleon(neutronMomentum, G4Neutron::Neutron() );
        G4InuclElementaryParticle spectatorNucleon(protonMomentum, G4Proton::Proton() );

        ScatteringProducts products = SingleNucleonScattering(projectile, targetNucleon);

        // Particles from SingleNucleonScattering are in CM frame of projectile
        // and moving neutron.  Transform back to lab frame with -betacm, then add
        // them to outgoing list
        globalOutput.reset();
        G4LorentzVector temp;
        for (G4int i = 0; i < G4int(products.size()); i++) {
          temp = products[i].getMomentum();
          temp.boost(-betacm);
          products[i].setMomentum(temp);
          globalOutput.addOutgoingParticle(products[i]);
        }

        // Add the recoil nucleon unmodified
        globalOutput.addOutgoingParticle(spectatorNucleon);

      } else {
        NucleonPair products = AbsorptionOnDeuteron(bullet);
        globalOutput.reset();
        globalOutput.addOutgoingParticle(products.first);
        globalOutput.addOutgoingParticle(products.second);
      }
    }  // Energy above threshold ?

    // Test code
    // G4int numPart = globalOutput.numberOfOutgoingParticles();
    // std::vector<G4InuclElementaryParticle> testList = globalOutput.getOutgoingParticles();
    // G4LorentzVector sumP;
    // G4cout << " Global output " << G4endl;
    // for (G4int i = 0; i < numPart; i++) {
    //  sumP += testList[i].getMomentum();
    //  G4cout << testList[i] << G4endl;
    // }
    // G4cout << " Global 4-momentum sum = " << sumP << G4endl;
    // G4cout << " Initial lab energy = " << mD + bullet->getEnergy() << G4endl;

  } else {
    G4Exception("G4LightTargetCollider::collide()","HAD_BERT_203",
                FatalException, "Scattering from this target not implemented");
  } 

  return;
}


G4double G4LightTargetCollider::GammaDCrossSection(G4double gammaEnergy)
{
  // Gamma deuteron cross section in mb parameterized from JLab data
  // No parameterization needed below pi0 threshold where cross section
  // is 100% disintegration
  G4double sigma = 1000.0;
  G4double term = 0.;
  if (gammaEnergy > 0.144 && gammaEnergy < 0.42) {
    term = (gammaEnergy - 0.24)/0.155;
    sigma = 0.065*std::exp(-term*term);
  } else if (gammaEnergy >= 0.42) {
    sigma = 0.000526/gammaEnergy/gammaEnergy/gammaEnergy/gammaEnergy;
  }

  return sigma;
}


NucleonPair G4LightTargetCollider::AbsorptionOnDeuteron(G4InuclParticle* bullet)
{
  // Do break-up in center of mass, convert to lab frame before returning 
  // particles

  G4double bulletMass = bullet->getMass();
  G4double bulletE = bullet->getEnergy();
 
  G4double S = bulletMass*bulletMass + mD*mD + 2.*mD*bulletE;

  G4double qcm = 0.;
  G4int outType1 = 0;
  G4int outType2 = 0;
  G4LorentzVector Mom1;
  G4LorentzVector Mom2;

  // Set up outgoing particle types
  if (bullet->getDefinition() == G4Gamma::Gamma() || 
      bullet->getDefinition() == G4PionZero::PionZero() ) {
    qcm = std::sqrt( (S - (mP + mN)*(mP + mN)) * (S - (mP - mN)*(mP - mN))/S/4.);
    Mom1.setE(std::sqrt(mP*mP + qcm*qcm) );
    outType1 = G4InuclParticleNames::proton;
    Mom2.setE(std::sqrt(mN*mN + qcm*qcm) );
    outType2 = G4InuclParticleNames::neutron;

  } else if (bullet->getDefinition() == G4PionPlus::PionPlus() ) {
    qcm = std::sqrt( (S - 4.*mP*mP)/4.);
    Mom1.setE(std::sqrt(mP*mP + qcm*qcm) );
    outType1 = G4InuclParticleNames::proton;
    Mom2.setE(std::sqrt(mP*mP + qcm*qcm) );
    outType2 = G4InuclParticleNames::proton;

  } else if (bullet->getDefinition() == G4PionMinus::PionMinus() ) {
    qcm = std::sqrt( (S - 4.*mN*mN)/4.);
    Mom1.setE(std::sqrt(mN*mN + qcm*qcm) );
    outType1 = G4InuclParticleNames::neutron;
    Mom2.setE(std::sqrt(mN*mN + qcm*qcm) );
    outType2 = G4InuclParticleNames::neutron;

  } else {
    G4Exception("G4LightTargetCollider::collide()","HAD_BERT_204",
                FatalException, "Illegal bullet type");
  } 

  // Sample angular distribution, assuming 100% S wave (no D-wave)
  G4ThreeVector qVect = qcm*G4RandomDirection();
  Mom1.setVect(qVect);
  Mom2.setVect(-qVect);

  // Boost to lab frame
  G4ThreeVector betacm(0., 0., bullet->getMomModule()/(bulletE + mD) );
  Mom1.boost(betacm);
  Mom2.boost(betacm);

  G4InuclElementaryParticle particle1(Mom1, outType1);
  G4InuclElementaryParticle particle2(Mom2, outType2);
  NucleonPair nucleon_pair(particle1, particle2); 

  // if pion, use parameterization of B.G. Ritchie, PRC 44, 533 (1991)
  // Total cross section: 1/E + Lorentzian

  return nucleon_pair;
}


ScatteringProducts
G4LightTargetCollider::SingleNucleonScattering(const G4InuclElementaryParticle& projectile,
                                               const G4InuclElementaryParticle& nucleon)
{
  // At this point projectile and nucleon momenta are in nucleon rest frame
  G4int reactionIndex = G4InuclElementaryParticle::type(projectile.getDefinition() ) 
                      * G4InuclElementaryParticle::type(nucleon.getDefinition() );

  const G4CascadeChannel* xsecTable = G4CascadeChannelTables::GetTable(reactionIndex);
  G4double ke = projectile.getKineticEnergy();
  G4int mult = xsecTable->getMultiplicity(ke);

  std::vector<G4double> masses;
  G4double mass = 0.0;
  G4LorentzVector totalMom = projectile.getMomentum() + nucleon.getMomentum();
  G4double Ecm = totalMom.mag(); 

  std::vector<G4LorentzVector> cmMomenta;
  std::vector<G4int> particle_kinds;
  G4int itry = 0;
  G4int itry_max = 200;
  G4bool generate = true;

  while (mult > 1) {
    itry = 0;
    generate = true;
    while (generate && itry < itry_max) {
      particle_kinds.clear();
      xsecTable->getOutgoingParticleTypes(particle_kinds, mult, ke);
      masses.clear();
      for (G4int i = 0; i < mult; i++) {
        mass = G4InuclElementaryParticle::getParticleMass(particle_kinds[i]);
        masses.push_back(mass);
      }

      fsGen.Configure(const_cast<G4InuclElementaryParticle*>(&projectile),
                      const_cast<G4InuclElementaryParticle*>(&nucleon),
                      particle_kinds);
      // Generate final state in CM of projectile and at-rest nucleon
      cmMomenta.clear();
      generate = !fsGen.Generate(Ecm, masses, cmMomenta);
      itry++;
    } // while 

    if (itry == itry_max) mult--;
    else break;

  } // while mult

  ScatteringProducts finalState;
  if (mult < 2) { 
    G4Exception("G4LightTargetCollider::SingleNucleonScattering()","HAD_BERT_202",
                 JustWarning, "Failed to generate final state");
    // Final state particles not in CM - just using them as dummies 
    finalState.push_back(projectile);
    finalState.push_back(nucleon);

  } else {
    for (G4int i = 0; i < mult; i++) {
      G4InuclElementaryParticle fsPart(cmMomenta[i], particle_kinds[i]);
      finalState.push_back(fsPart);
    }
  }

  return finalState;
}

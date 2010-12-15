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
#include "G4InclInput.hh"

G4InclInput::G4InclInput(const G4HadProjectile &aTrack, const G4Nucleus &theNucleus, G4bool inverseKinematics = false) {
  usingInverseKinematics = inverseKinematics;

  fTargetA = theNucleus.GetA_asInt(); // Target mass number
  fTargetZ = theNucleus.GetZ_asInt(); // Target charge number
  fBulletType = getBulletType(aTrack.GetDefinition()); // Projectile type (INCL particle code)
  fBulletE = aTrack.GetKineticEnergy() / MeV; // Projectile energy (total, in MeV)
  fTimeScale = 1.0; // Time scaling
  fNuclearPotential = 45.0; // Nuclear potential
  setExtendedProjectileInfo(aTrack.GetDefinition());
  icoup = 0;
  breakupThreshold = 10;
  fMinNeutronEnergy = 0.0;
  fMinProtonE = 0.0;
}

G4InclInput::~G4InclInput() {}

void G4InclInput::printInfo() {
  G4cout <<"Target: A = " << targetA() << " Z = " << targetZ() << G4endl;
  G4cout <<"Projectile: type = " << bulletType() << " energy = " << bulletE() << G4endl;
}

void G4InclInput::printProjectileTargetInfo(const G4HadProjectile &aTrack, const G4Nucleus &theNucleus) {
  G4cout <<"Projectile = " << aTrack.GetDefinition()->GetParticleName() << G4endl;
  G4cout <<"    four-momentum: " << aTrack.Get4Momentum() << G4endl;
  G4cout <<"Energy = " << aTrack.GetKineticEnergy() / MeV << G4endl;
  G4cout <<"Target A = " << theNucleus.GetA_asInt() << " Z = " << theNucleus.GetZ_asInt() << G4endl;
}


G4bool G4InclInput::canUseInverseKinematics(const G4HadProjectile &aTrack, const G4Nucleus &theNucleus) {
  G4int targetA = theNucleus.GetA_asInt();
  const G4ParticleDefinition *projectileDef = aTrack.GetDefinition();
  G4int projectileA = projectileDef->GetAtomicMass();
  //    G4int projectileZ = projectileDef->GetAtomicNumber();
  if(targetA > 0 && targetA < 18 && (projectileDef != G4Proton::Proton() &&
				     projectileDef != G4Neutron::Neutron() &&
				     projectileDef != G4PionPlus::PionPlus() &&
				     projectileDef != G4PionZero::PionZero() &&
				     projectileDef != G4PionMinus::PionMinus()) &&
			      projectileA > 1) {
    return true;
  } else {
    return false;
  }
}

void G4InclInput::setExtendedProjectileInfo(const G4ParticleDefinition *pd) {
  if(getBulletType(pd) == -666) {
    theExtendedProjectileA = pd->GetAtomicMass();
    theExtendedProjectileZ = pd->GetAtomicNumber();
    isExtended = true;
  } else {
    isExtended = false;
  }
}

G4int G4InclInput::getBulletType(const G4ParticleDefinition *pd) {
  //  G4ParticleTable *pt = G4ParticleTable::GetParticleTable();

  if(pd == G4Proton::Proton()) {
    return 1;
  } else if(pd == G4Neutron::Neutron()) {
    return 2;
  } else if(pd == G4PionPlus::PionPlus()) {
    return 3;
  } else if(pd == G4PionMinus::PionMinus()) {
    return 5;
  } else if(pd == G4PionZero::PionZero()) {
    return 4;
  } else if(pd == G4Deuteron::Deuteron()) {
    return 6;
  } else if(pd == G4Triton::Triton()) {
    return 7;
  } else if(pd == G4He3::He3()) {
    return 8;
  } else if(pd == G4Alpha::Alpha()) {
    return 9;
    //  } else if(pd == pt->GetIon(6, 12, 0.0)) { // C12 special case. This should be phased-out in favor of "extended projectile"
    //    return -12;
  } else { // Is this extended projectile?
    G4int A = pd->GetAtomicMass();
    G4int Z = pd->GetAtomicNumber();
    if(A > 4 && A <= 16 && Z > 2 && Z <= 8) { // Ions from Lithium to Oxygen
      return -666; // Code of an extended projectile
    }
  }
  G4cout <<"Error! Projectile " << pd->GetParticleName() << " not defined!" << G4endl;
  return 0;
}

G4ParticleDefinition* G4InclInput::getParticleDefinition(G4int inclParticleCode) {
  switch(inclParticleCode) {
  case 1:
    return G4Proton::ProtonDefinition();
    break;
  case 2:
    return G4Neutron::NeutronDefinition();
    break;
  case 3:
    return G4PionPlus::PionPlusDefinition();
    break;
  case 4:
    return G4PionMinus::PionMinusDefinition();
    break;
  case 5:
    return G4PionZero::PionZeroDefinition();
    break;
  case 6:
    return G4Deuteron::DeuteronDefinition();
    break;
  case 7:
    return G4Triton::Triton();
    break;
  case 8:
    return G4He3::He3Definition();
    break;
  case 9:
    return G4Alpha::AlphaDefinition();
    break;
  }
  return 0;
}

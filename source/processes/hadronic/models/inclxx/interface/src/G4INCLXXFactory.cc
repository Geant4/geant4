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
// INCL++ revision: v5.0_rc3
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4INCLXXFactory.hh"
#include "G4ParticleTable.hh"

G4INCL::ParticleType G4INCLXXFactory::toINCLParticleType(const G4ParticleDefinition *pdef) {
  if(     pdef == G4Proton::Proton())       return G4INCL::Proton;
  else if(pdef == G4Neutron::Neutron())     return G4INCL::Neutron;
  else if(pdef == G4PionPlus::PionPlus())   return G4INCL::PiPlus;
  else if(pdef == G4PionMinus::PionMinus()) return G4INCL::PiMinus;
  else if(pdef == G4PionZero::PionZero())   return G4INCL::PiZero;
  else                                      return G4INCL::UnknownParticle;
}

const G4ParticleDefinition* G4INCLXXFactory::fromINCLParticleType(G4INCL::ParticleType ptype) {
  if(     ptype == G4INCL::Proton)          return G4Proton::Proton();
  else if(ptype == G4INCL::Neutron)         return G4Neutron::Neutron();
  else if(ptype == G4INCL::PiPlus)          return G4PionPlus::PionPlus();
  else if(ptype == G4INCL::PiMinus)         return G4PionMinus::PionMinus();
  else if(ptype == G4INCL::PiZero)          return G4PionZero::PionZero();
  else if(ptype == G4INCL::UnknownParticle) return 0;
  else                                      return 0;
}

G4INCL::Particle* G4INCLXXFactory::createProjectile(const G4HadProjectile &aTrack) {
  const G4ParticleDefinition *pdef = aTrack.GetDefinition();
  G4INCL::ParticleType projectileType = G4INCLXXFactory::toINCLParticleType(pdef);
  const G4double kineticEnergy = aTrack.GetKineticEnergy();
  const G4double mass = G4INCL::ParticleTable::getMass(projectileType);
  const G4double energy = kineticEnergy + mass;
  const G4double pz = std::sqrt(energy*energy - mass*mass);
  G4INCL::ThreeVector momentum(0.0, 0.0, pz);
  G4INCL::ThreeVector position(0.0, 0.0, 0.0); // Projectile position
					       // doesn't actually
					       // matter.

  G4INCL::Particle *projectile = new G4INCL::Particle(projectileType, energy,
						      momentum, position);
  return projectile;
}

G4INCL::INCL* G4INCLXXFactory::createModel(const G4Nucleus &theNucleus) {
  G4int A = theNucleus.GetA_asInt();
  G4int Z = theNucleus.GetZ_asInt();
  G4INCL::Config *theConfig = new G4INCL::Config(A, Z, G4INCL::Proton, 1200.0);
  theConfig->setTargetA(A);
  theConfig->setTargetZ(Z);

  G4INCL::INCL *theINCLModel = new G4INCL::INCL(theConfig);

  return theINCLModel;
}

G4ParticleDefinition* G4INCLXXFactory::toG4ParticleDefinition(G4int A,
							      G4int Z) {
  if     (A == 1 && Z == 1)  return G4Proton::Proton();
  else if(A == 1 && Z == 0)  return G4Neutron::Neutron();
  else if(A == 0 && Z == 1)  return G4PionPlus::PionPlus();
  else if(A == 0 && Z == -1) return G4PionMinus::PionMinus();
  else if(A == 0 && Z == 0)  return G4PionZero::PionZero();
  else if(A == 2 && Z == 1)  return G4Deuteron::Deuteron();
  else if(A == 3 && Z == 1)  return G4Triton::Triton();
  else if(A == 3 && Z == 2)  return G4He3::He3();
  else if(A == 4 && Z == 2)  return G4Alpha::Alpha();
  else if(A > 0 && Z > 0 && A > Z) { // Returns ground state ion definition
    return G4ParticleTable::GetParticleTable()->GetIon(Z, A, 0.0);
  } else { // Error, unrecognized particle
    return 0;
  }
}

G4DynamicParticle* G4INCLXXFactory::toG4Particle(G4int A, G4int Z,
						 G4double kinE,
						 G4double px,
						 G4double py, G4double pz) {
  const G4ParticleDefinition *def = toG4ParticleDefinition(A, Z);
  if(def == 0) { // Check if we have a valid particle definition
    return 0;
  }
  const G4double energy = kinE / MeV;
  const G4ThreeVector momentum(px, py, pz);
  const G4ThreeVector momentumDirection = momentum.unit();
  G4DynamicParticle *p = new G4DynamicParticle(def, momentumDirection, energy);
  return p;
}

G4double G4INCLXXFactory::remnant4MomentumScaling(G4double mass,
						  G4double kineticE,
						  G4double px, G4double py,
						  G4double pz) {
  const G4double p2 = px*px + py*py + pz*pz;
  if(p2 > 0.0) {
    const G4double pnew2 = kineticE*kineticE + 2.0*kineticE*mass;
    return std::sqrt(pnew2)/std::sqrt(p2);
  } else {
    return 1.0;
  }
}


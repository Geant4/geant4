//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4CollisionNNToNDelta1905.cc,v 1.1 2003/10/07 12:37:36 hpw Exp $ //

#include "globals.hh"
#include "G4CollisionNNToNDelta1905.hh"
#include "G4ConcreteNNToNDeltaStar.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4ParticleTable.hh"
#include "G4ShortLivedConstructor.hh"

// complete hpw

G4CollisionNNToNDelta1905::G4CollisionNNToNDelta1905()
{ 
  // Subtype of interacting particles
  G4ShortLivedConstructor ShortLived;
  ShortLived.ConstructParticle();

  G4ParticleDefinition * aProton = G4Proton::ProtonDefinition();
  G4ParticleDefinition * aNeutron = G4Neutron::NeutronDefinition();
  
  G4ParticleDefinition * aDm_1905 = G4ParticleTable::GetParticleTable()->FindParticle(1116); 
  G4ParticleDefinition * aD0_1905 = G4ParticleTable::GetParticleTable()->FindParticle(1216); 
  G4ParticleDefinition * aDp_1905 = G4ParticleTable::GetParticleTable()->FindParticle(2126); 
  G4ParticleDefinition * aDpp_1905 = G4ParticleTable::GetParticleTable()->FindParticle(2226); 
  G4CollisionComposite::AddComponent(new G4ConcreteNNToNDeltaStar(aNeutron, aNeutron, aNeutron, aD0_1905));  
  G4CollisionComposite::AddComponent(new G4ConcreteNNToNDeltaStar(aNeutron, aNeutron, aProton, aDm_1905));
  G4CollisionComposite::AddComponent(new G4ConcreteNNToNDeltaStar(aNeutron, aProton, aProton, aD0_1905));
  G4CollisionComposite::AddComponent(new G4ConcreteNNToNDeltaStar(aNeutron, aProton, aNeutron, aDp_1905));
  G4CollisionComposite::AddComponent(new G4ConcreteNNToNDeltaStar(aProton, aProton, aNeutron, aDpp_1905));
  G4CollisionComposite::AddComponent(new G4ConcreteNNToNDeltaStar(aProton, aProton, aProton, aDp_1905));
}




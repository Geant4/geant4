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
// $Id: G4CollisionNNToDeltaDelta.cc,v 1.1 2003/10/07 12:37:34 hpw Exp $ //

#include "globals.hh"
#include "G4CollisionNNToDeltaDelta.hh"
#include "G4KineticTrack.hh"
#include "G4VCrossSectionSource.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4XAqmElastic.hh"
#include "G4AngularDistribution.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
#include "G4KineticTrackVector.hh"
#include "G4ParticleTable.hh"
#include "G4CollisionVector.hh"
#include "G4ConcreteNNToDeltaDelta.hh"

G4CollisionNNToDeltaDelta::G4CollisionNNToDeltaDelta()
{ 
  G4ParticleDefinition * aProton = G4Proton::ProtonDefinition();
  G4ParticleDefinition * aNeutron = G4Neutron::NeutronDefinition();
  G4ParticleDefinition * aDeltapp = G4ParticleTable::GetParticleTable()->FindParticle(2224); 
  G4ParticleDefinition * aDeltap = G4ParticleTable::GetParticleTable()->FindParticle(2214); 
  G4ParticleDefinition * aDelta0 = G4ParticleTable::GetParticleTable()->FindParticle(2114); 
  G4ParticleDefinition * aDeltam = G4ParticleTable::GetParticleTable()->FindParticle(1114); 
  
  G4CollisionComposite::AddComponent(new G4ConcreteNNToDeltaDelta(aNeutron, aNeutron, aDelta0, aDelta0));  
  G4CollisionComposite::AddComponent(new G4ConcreteNNToDeltaDelta(aNeutron, aNeutron, aDeltam, aDeltap));  
  G4CollisionComposite::AddComponent(new G4ConcreteNNToDeltaDelta(aNeutron, aProton, aDelta0, aDeltap));  
  G4CollisionComposite::AddComponent(new G4ConcreteNNToDeltaDelta(aNeutron, aProton, aDeltam, aDeltapp));  
  G4CollisionComposite::AddComponent(new G4ConcreteNNToDeltaDelta(aProton,  aProton, aDeltap, aDeltap));  
  G4CollisionComposite::AddComponent(new G4ConcreteNNToDeltaDelta(aProton,  aProton, aDelta0, aDeltapp));  

}


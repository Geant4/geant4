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
// $Id: G4CollisionNNToNDelta.cc,v 1.3 2003-12-12 12:28:08 hpw Exp $ //

#include "globals.hh"
#include "G4CollisionNNToNDelta.hh"
#include "G4ConcreteNNToNDelta.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4ParticleTable.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4Pair.hh"
#include "G4HadParticleCodes.hh"

// complete hpw
typedef INT4(G4ConcreteNNToNDelta, NeutronPC, NeutronPC, NeutronPC, Delta0PC)  theC1;
typedef INT4(G4ConcreteNNToNDelta, NeutronPC, NeutronPC, ProtonPC,  DeltamPC)  theC2;
typedef INT4(G4ConcreteNNToNDelta, NeutronPC, ProtonPC,  ProtonPC,  Delta0PC)  theC3;
typedef INT4(G4ConcreteNNToNDelta, NeutronPC, ProtonPC,  NeutronPC, DeltapPC)  theC4;
typedef INT4(G4ConcreteNNToNDelta, ProtonPC,  ProtonPC,  NeutronPC, DeltappPC) theC5;
typedef INT4(G4ConcreteNNToNDelta, ProtonPC,  ProtonPC,  ProtonPC,  DeltapPC)  theC6;

typedef GROUP6(theC1, theC2, theC3, theC4, theC5, theC6) theChannels;
       
G4CollisionNNToNDelta::G4CollisionNNToNDelta()
{ 
  // Subtype of interacting particles
  G4ShortLivedConstructor ShortLived;
  ShortLived.ConstructParticle();
  G4ForEach<theChannels, Resolve>::Apply(this);

}




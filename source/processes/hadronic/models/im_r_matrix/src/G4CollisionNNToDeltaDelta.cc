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
// $Id: G4CollisionNNToDeltaDelta.cc,v 1.2.2.1 2004/03/24 13:18:30 hpw Exp $ //

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
#include "G4Pair.hh"
#include "G4HadParticleCodes.hh"

typedef G4ConcreteNNToDeltaDelta channelType;
typedef INT4<channelType, NeutronPC, NeutronPC, Delta0PC, Delta0PC>  theC1;
typedef INT4<channelType, NeutronPC, NeutronPC, DeltamPC, DeltapPC>  theC2;
typedef INT4<channelType, NeutronPC, ProtonPC,  Delta0PC, DeltapPC>  theC3;
typedef INT4<channelType, NeutronPC, ProtonPC,  DeltamPC, DeltappPC> theC4;
typedef INT4<channelType, ProtonPC,  ProtonPC,  DeltapPC, DeltapPC>  theC5;
typedef INT4<channelType, ProtonPC,  ProtonPC,  Delta0PC, DeltappPC> theC6;

typedef GROUP6(theC1, theC2, theC3, theC4, theC5, theC6) theChannels;

G4CollisionNNToDeltaDelta::G4CollisionNNToDeltaDelta()
{ 
  Resolve aR;
  G4ForEach<theChannels>::Apply(&aR, this);
}


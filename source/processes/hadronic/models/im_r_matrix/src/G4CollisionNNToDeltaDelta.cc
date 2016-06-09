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
//
// $Id: G4CollisionNNToDeltaDelta.cc,v 1.3 2010-03-12 15:45:18 gunter Exp $ //

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


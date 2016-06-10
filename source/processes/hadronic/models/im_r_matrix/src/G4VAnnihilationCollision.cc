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

#include "globals.hh"
#include "G4VAnnihilationCollision.hh"
#include "G4KineticTrack.hh"
#include "G4VCrossSectionSource.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4XNNElastic.hh"
#include "G4AngularDistribution.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
#include "G4KineticTrackVector.hh"
#include "Randomize.hh"
#include "G4PionPlus.hh"


G4VAnnihilationCollision::G4VAnnihilationCollision()
{ 
	theAngularDistribution=0;
}


G4VAnnihilationCollision::~G4VAnnihilationCollision()
{
	if (theAngularDistribution) {
		delete theAngularDistribution;
		theAngularDistribution=0;
	}
}


G4KineticTrackVector* G4VAnnihilationCollision::FinalState(const G4KineticTrack& trk1, 
							      const G4KineticTrack& trk2) const
{ 
  G4LorentzVector p = trk1.Get4Momentum() + trk2.Get4Momentum();
  G4double sqrtS = p.m();
//  G4double s = sqrtS * sqrtS;

//  G4double m1 = trk1.GetActualMass();
//  G4double m2 = trk2.GetActualMass();


  // Unit vector of three-momentum

  G4ThreeVector pFinal(0.0, 0.0, 0.0);
  G4double eFinal = sqrtS;
 
  G4LorentzVector p4Final(pFinal, eFinal); 

  // Lorentz transformation
  G4LorentzRotation toLabFrame(p.boostVector());
  p4Final *= toLabFrame;
  
 const G4ParticleDefinition* OutputDefinition = GetOutgoingParticle(trk1,trk2);
 G4KineticTrack* final = new G4KineticTrack(OutputDefinition, 0.0, trk1.GetPosition(), p4Final);
  
  G4KineticTrackVector* finalTracks = new G4KineticTrackVector;

  finalTracks->push_back(final);

  return finalTracks;
}


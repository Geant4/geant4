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
}


G4VAnnihilationCollision::~G4VAnnihilationCollision()
{ 
}


G4KineticTrackVector* G4VAnnihilationCollision::FinalState(const G4KineticTrack& trk1, 
							      const G4KineticTrack& trk2) const
{ 
  G4LorentzVector p = trk1.Get4Momentum() + trk2.Get4Momentum();
  G4double sqrtS = p.m();
//  G4double s = sqrtS * sqrtS;

//  G4double m1 = trk1.GetActualMass();
//  G4double m2 = trk2.GetActualMass();

  G4ParticleDefinition* OutputDefinition = const_cast<G4ParticleDefinition *>( GetOutgoingParticle(trk1,trk2) );

  // Unit vector of three-momentum

  G4ThreeVector pFinal(0.0, 0.0, 0.0);
  G4double eFinal = sqrtS;
 
  G4LorentzVector p4Final(pFinal, eFinal); 

  // Lorentz transformation
  G4LorentzRotation toLabFrame(p.boostVector());
  p4Final *= toLabFrame;
  
  G4KineticTrack* final = new G4KineticTrack(OutputDefinition, 0.0, 
					     trk1.GetPosition(), p4Final);
  
  G4KineticTrackVector* finalTracks = new G4KineticTrackVector;

  finalTracks->push_back(final);

  return finalTracks;
}


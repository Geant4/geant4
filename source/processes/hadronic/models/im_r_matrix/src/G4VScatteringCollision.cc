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
// @hpw@ misses the sampling of two breit wigner in a corelated fashion, 
// @hpw@ to be usefull for resonance resonance scattering.

#include <typeinfo>

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4VScatteringCollision.hh"
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

G4VScatteringCollision::G4VScatteringCollision()
{ 
  theAngularDistribution = new G4AngularDistribution(true);
}


G4VScatteringCollision::~G4VScatteringCollision()
{ 
  delete theAngularDistribution;
  theAngularDistribution=0;
}


G4KineticTrackVector* G4VScatteringCollision::FinalState(const G4KineticTrack& trk1, 
							    const G4KineticTrack& trk2) const
{ 
  const G4VAngularDistribution* angDistribution = GetAngularDistribution();
  G4LorentzVector p = trk1.Get4Momentum() + trk2.Get4Momentum();
  G4double sqrtS = p.m();
  G4double S = sqrtS * sqrtS;

  std::vector<const G4ParticleDefinition*> OutputDefinitions = GetOutgoingParticles();
  if (OutputDefinitions.size() != 2)
    throw G4HadronicException(__FILE__, __LINE__, "G4VScatteringCollision: Too many output particles!");

  if (OutputDefinitions[0]->IsShortLived() && OutputDefinitions[1]->IsShortLived())
  {
    if(getenv("G4KCDEBUG")) G4cerr << "two shortlived for Type = "<<typeid(*this).name()<<G4endl;
    // throw G4HadronicException(__FILE__, __LINE__, "G4VScatteringCollision: can't handle two shortlived particles!"); // @hpw@
  }
  
  G4double outm1 = OutputDefinitions[0]->GetPDGMass();
  G4double outm2 = OutputDefinitions[1]->GetPDGMass();

  if (OutputDefinitions[0]->IsShortLived())
  {
    outm1 = SampleResonanceMass(outm1, 
                OutputDefinitions[0]->GetPDGWidth(),
		G4Neutron::NeutronDefinition()->GetPDGMass()+G4PionPlus::PionPlus()->GetPDGMass(),
		sqrtS-(G4Neutron::NeutronDefinition()->GetPDGMass()+G4PionPlus::PionPlus()->GetPDGMass()));

  }
  if (OutputDefinitions[1]->IsShortLived())
  {
    outm2 = SampleResonanceMass(outm2, OutputDefinitions[1]->GetPDGWidth(),
			G4Neutron::NeutronDefinition()->GetPDGMass()+G4PionPlus::PionPlus()->GetPDGMass(),
			sqrtS-outm1);
  }
  
  // Angles of outgoing particles
  G4double cosTheta = angDistribution->CosTheta(S, trk1.GetActualMass(), trk2.GetActualMass());
  G4double phi = angDistribution->Phi();

  // Unit vector of three-momentum
  G4LorentzRotation fromCMSFrame(p.boostVector());
  G4LorentzRotation toCMSFrame(fromCMSFrame.inverse());
  G4LorentzVector TempPtr = toCMSFrame*trk1.Get4Momentum();
  G4LorentzRotation toZ;
  toZ.rotateZ(-1*TempPtr.phi());
  toZ.rotateY(-1*TempPtr.theta());
  G4LorentzRotation toCMS(toZ.inverse());

  G4ThreeVector pFinal1(std::sin(std::acos(cosTheta))*std::cos(phi), std::sin(std::acos(cosTheta))*std::sin(phi), cosTheta);

  // Three momentum in cm system
  G4double pCM = std::sqrt( (S-(outm1+outm2)*(outm1+outm2)) * (S-(outm1-outm2)*(outm1-outm2)) /(4.*S));
  pFinal1 = pFinal1 * pCM;
  G4ThreeVector pFinal2 = -pFinal1;

  G4double eFinal1 = std::sqrt(pFinal1.mag2() + outm1*outm1);
  G4double eFinal2 = std::sqrt(pFinal2.mag2() + outm2*outm2);

  G4LorentzVector p4Final1(pFinal1, eFinal1);
  G4LorentzVector p4Final2(pFinal2, eFinal2);
  p4Final1 = toCMS*p4Final1;
  p4Final2 = toCMS*p4Final2;


  // Lorentz transformation
  G4LorentzRotation toLabFrame(p.boostVector());
  p4Final1 *= toLabFrame;
  p4Final2 *= toLabFrame;

  // Final tracks are copies of incoming ones, with modified 4-momenta

  G4double chargeBalance = OutputDefinitions[0]->GetPDGCharge()+OutputDefinitions[1]->GetPDGCharge();
  chargeBalance-= trk1.GetDefinition()->GetPDGCharge();
  chargeBalance-= trk2.GetDefinition()->GetPDGCharge();
  if(std::abs(chargeBalance) >.1)
  {
    G4cout << "Charges in "<<typeid(*this).name()<<G4endl;
    G4cout << OutputDefinitions[0]->GetPDGCharge()<<" "<<OutputDefinitions[0]->GetParticleName()
           << OutputDefinitions[1]->GetPDGCharge()<<" "<<OutputDefinitions[1]->GetParticleName()
	   << trk1.GetDefinition()->GetPDGCharge()<<" "<<trk1.GetDefinition()->GetParticleName()
	   << trk2.GetDefinition()->GetPDGCharge()<<" "<<trk2.GetDefinition()->GetParticleName()<<G4endl;
  }
  G4KineticTrack* final1 = new G4KineticTrack(OutputDefinitions[0], 0.0, trk1.GetPosition(), p4Final1);
  G4KineticTrack* final2 = new G4KineticTrack(OutputDefinitions[1], 0.0, trk2.GetPosition(), p4Final2);

  G4KineticTrackVector* finalTracks = new G4KineticTrackVector;

  finalTracks->push_back(final1);
  finalTracks->push_back(final2);

  return finalTracks;
}



double G4VScatteringCollision::SampleResonanceMass(const double poleMass, 
						   const double gamma,
						   const double aMinMass,
						   const double maxMass) const
{
  // Chooses a mass randomly between minMass and maxMass 
  //     according to a Breit-Wigner function with constant 
  //     width gamma and pole poleMass

  G4double minMass = aMinMass;
  if (minMass > maxMass) G4cerr << "##################### SampleResonanceMass: particle out of mass range" << G4endl;
  if(minMass > maxMass) minMass -= G4PionPlus::PionPlus()->GetPDGMass();
  if(minMass > maxMass) minMass = 0;

  if (gamma < 1E-10*GeV)
    return std::max(minMass,std::min(maxMass, poleMass));
  else {
    double fmin = BrWigInt0(minMass, gamma, poleMass);
    double fmax = BrWigInt0(maxMass, gamma, poleMass);
    double f = fmin + (fmax-fmin)*G4UniformRand();
    return BrWigInv(f, gamma, poleMass);
  }
}

void G4VScatteringCollision::establish_G4MT_TLS_G4VScatteringCollision()
{
  establish_G4MT_TLS_G4VCollision();
  if ( theAngularDistribution ) delete theAngularDistribution;
  theAngularDistribution = new G4AngularDistribution(true);
}

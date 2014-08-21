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
// $Id: G4ConcreteNNTwoBodyResonance.cc,v 1.2 2006-06-29 20:40:00 gunter Exp $ //

#include "globals.hh"
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
#include "G4XResonance.hh"
#include "G4ConcreteNNTwoBodyResonance.hh"

G4ConcreteNNTwoBodyResonance::G4ConcreteNNTwoBodyResonance(const G4ParticleDefinition* aPrimary,
					   const G4ParticleDefinition* bPrimary,
					   const G4ParticleDefinition* aSecondary,
					   const G4ParticleDefinition* bSecondary,
		                           const G4VXResonanceTable& sigmaTable)
  : thePrimary1(aPrimary), thePrimary2(bPrimary)
{
  theOutGoing.push_back(aSecondary);
  theOutGoing.push_back(bSecondary);

  crossSectionSource = new G4XResonance(aPrimary, bPrimary, 
					aSecondary->GetPDGiIsospin(), 
					aSecondary->GetPDGiSpin(),
					aSecondary->GetPDGMass(),
					bSecondary->GetPDGiIsospin(), 
					bSecondary->GetPDGiSpin(),
					bSecondary->GetPDGMass(),
					aSecondary->GetParticleName(),
					bSecondary->GetParticleName(),
					sigmaTable);
}

G4ConcreteNNTwoBodyResonance::~G4ConcreteNNTwoBodyResonance()
{ 
  if (crossSectionSource) delete crossSectionSource;
  crossSectionSource=0;
}

G4bool G4ConcreteNNTwoBodyResonance::IsInCharge(const G4KineticTrack& trk1, 
					 const G4KineticTrack& trk2) const
{
  if (trk1.GetDefinition()==thePrimary1 && trk2.GetDefinition()==thePrimary2) return true;
  if (trk1.GetDefinition()==thePrimary2 && trk2.GetDefinition()==thePrimary1) return true;
  return false;
}

G4ConcreteNNTwoBodyResonance::G4ConcreteNNTwoBodyResonance(void *, void *, void *, void *, void *, void *, void *)
    : crossSectionSource(0), thePrimary1(0), thePrimary2(0)
{}

void G4ConcreteNNTwoBodyResonance::establish_G4MT_TLS_G4ConcreteNNTwoBodyResonance(const G4ParticleDefinition* aPrimary,
					   const G4ParticleDefinition* bPrimary,
					   const G4ParticleDefinition* aSecondary,
					   const G4ParticleDefinition* bSecondary,
		               const G4VXResonanceTable& sigmaTable)
{
  establish_G4MT_TLS_G4VScatteringCollision();
  thePrimary1=aPrimary;
  thePrimary2=bPrimary;
  theOutGoing.push_back(aSecondary);
  theOutGoing.push_back(bSecondary);

  crossSectionSource = new G4XResonance(aPrimary, bPrimary, 
					aSecondary->GetPDGiIsospin(), 
					aSecondary->GetPDGiSpin(),
					aSecondary->GetPDGMass(),
					bSecondary->GetPDGiIsospin(), 
					bSecondary->GetPDGiSpin(),
					bSecondary->GetPDGMass(),
					aSecondary->GetParticleName(),
					bSecondary->GetParticleName(),
					sigmaTable);
}

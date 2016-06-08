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
// $Id: G4ConcreteNNTwoBodyResonance.cc,v 1.2 2002/12/12 19:17:50 gunter Exp $ //

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
  delete crossSectionSource;
}

G4bool G4ConcreteNNTwoBodyResonance::IsInCharge(const G4KineticTrack& trk1, 
					 const G4KineticTrack& trk2) const
{
  if (trk1.GetDefinition()==thePrimary1 && trk2.GetDefinition()==thePrimary2) return true;
  if (trk1.GetDefinition()==thePrimary2 && trk2.GetDefinition()==thePrimary1) return true;
  return false;
}


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
// $Id: G4HadronElastic.hh 90228 2015-05-21 08:49:57Z gcosmo $
//
// Geant4 Header : G4HadronElastic
//
// Author : V.Ivanchenko 29 June 2009 (redesign old elastic model)
//  
// Modified:
//
// Class Description
// Default model for elastic scattering; GHEISHA algorithm is used 
// Class Description - End

#ifndef G4HadronElastic_h
#define G4HadronElastic_h 1
 
#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4NucleiProperties.hh"

class G4ParticleDefinition;

class G4HadronElastic : public G4HadronicInteraction
{
public:

  G4HadronElastic(const G4String& name = "hElasticLHEP");

  virtual ~G4HadronElastic();
 
  // implementation of the G4HadronicInteraction interface
  virtual G4HadFinalState * ApplyYourself(const G4HadProjectile & aTrack, 
					  G4Nucleus & targetNucleus);

  // sample momentum transfer using Lab. momentum
  virtual G4double SampleInvariantT(const G4ParticleDefinition* p, 
				    G4double plab,
				    G4int Z, G4int A);

  inline void SetLowestEnergyLimit(G4double value);

  inline G4double LowestEnergyLimit() const;

  inline G4double ComputeMomentumCMS(const G4ParticleDefinition* p, 
				     G4double plab, G4int Z, G4int A);
  
  virtual void ModelDescription(std::ostream&) const;

private:

  G4ParticleDefinition* theProton;
  G4ParticleDefinition* theNeutron;
  G4ParticleDefinition* theDeuteron;
  G4ParticleDefinition* theAlpha;

  G4double lowestEnergyLimit;  

};

inline void G4HadronElastic::SetLowestEnergyLimit(G4double value)
{
  lowestEnergyLimit = value;
}

inline G4double G4HadronElastic::LowestEnergyLimit() const
{
  return lowestEnergyLimit;
}

inline G4double
G4HadronElastic::ComputeMomentumCMS(const G4ParticleDefinition* p, 
				     G4double plab, G4int Z, G4int A)
{
  G4double m1 = p->GetPDGMass();
  G4double m12= m1*m1;
  G4double mass2 = G4NucleiProperties::GetNuclearMass(A, Z);
  return plab*mass2/std::sqrt(m12 + mass2*mass2 + 2.*mass2*std::sqrt(m12 + plab*plab));
}

#endif

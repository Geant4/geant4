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
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4SampleResonance ----------------
//             by Henning Weber, March 2001.
//      helper class for sampling resonance masses
// ------------------------------------------------------------


#include "globals.hh"
#include <iostream>
#include "G4SampleResonance.hh"
#include "G4DecayTable.hh"
#include "Randomize.hh"
#include "G4HadronicException.hh"

G4SampleResonance::minMassMapType G4SampleResonance::minMassCache;

G4double G4SampleResonance::GetMinimumMass(const G4ParticleDefinition* p) const
{
  
  G4double minResonanceMass = DBL_MAX;

  if ( p->IsShortLived() ) 
    {
      minMassMapIterator i = minMassCache.find(p);
      if ( i!=minMassCache.end() ) 
	{
	  minResonanceMass = (*i).second; 
	} 
      else 
	{
	  // G4cout << "--> request for " << p->GetParticleName() << G4endl;

	  const G4DecayTable* theDecays = const_cast<G4ParticleDefinition*>(p)->GetDecayTable();
	  const G4int nDecays = theDecays->entries();
      
	  for (G4int i=0; i<nDecays; i++) 
	    {
	      const G4VDecayChannel* aDecay = theDecays->GetDecayChannel(i);
	      const G4int nDaughters = aDecay->GetNumberOfDaughters();
	      
	      G4double minChannelMass = 0;
	      
	      for (G4int j=0; j<nDaughters; j++) 
		{
		  const G4ParticleDefinition* aDaughter = const_cast<G4VDecayChannel*>(aDecay)->GetDaughter(j);
		  G4double minMass = GetMinimumMass(aDaughter);
		  if (!minMass) minMass = DBL_MAX; // exclude gamma channel;
		  minChannelMass+=minMass;
		}
	      // G4cout << "channel mass for the above is " << minChannelMass/MeV << G4endl;
	      if (minChannelMass < minResonanceMass) minResonanceMass = minChannelMass;
	      
	    }
	  // replace this as soon as the compiler supports mutable!!
	  G4SampleResonance* self = const_cast<G4SampleResonance*>(this);
	  (self->minMassCache)[p] = minResonanceMass;
	  
	}
    } 
  else 
    {
    
      minResonanceMass = p->GetPDGMass();
      
    }
  // G4cout << "minimal mass for " << p->GetParticleName() << " is " << minResonanceMass/MeV << G4endl;
  
  return minResonanceMass;
} 
      


G4double G4SampleResonance::SampleMass(const G4ParticleDefinition* p, const G4double maxMass) const
{
  return SampleMass(p->GetPDGMass(), p->GetPDGWidth(), GetMinimumMass(p), maxMass);
}


G4double G4SampleResonance::SampleMass(const G4double poleMass, 
				       const G4double gamma,
				       const G4double minMass,
				       const G4double maxMass) const
{
  // Chooses a mass randomly between minMass and maxMass 
  //     according to a Breit-Wigner function with constant 
  //     width gamma and pole poleMass

  if ( minMass > maxMass ) 
    {
      throw G4HadronicException(__FILE__, __LINE__,
           "SampleResonanceMass: mass range negative (minMass>maxMass)");
    }

  G4double returnMass;

  if ( gamma < DBL_EPSILON ) 
    {
      returnMass = std::max(minMass, std::min(maxMass, poleMass));
    } 
  else 
    {
      double fmin = BrWigInt0(minMass, gamma, poleMass);
      double fmax = BrWigInt0(maxMass, gamma, poleMass);
      double f = fmin + (fmax-fmin)*G4UniformRand();
      returnMass = BrWigInv(f, gamma, poleMass);      
    }

  return returnMass;
}





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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4SampleResonance.cc,v 1.4 2001/10/06 08:17:02 hpw Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4SampleResonance ----------------
//             by Henning Weber, March 2001.
//      helper class for sampling resonance masses
// ------------------------------------------------------------


#include "globals.hh"
#include "g4std/iostream"
#include "G4SampleResonance.hh"
#include "G4DecayTable.hh"
#include "Randomize.hh"


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
      G4Exception("SampleResonanceMass: mass range negative (minMass>maxMass)");
    }

  G4double returnMass;

  if ( gamma < DBL_EPSILON ) 
    {
      returnMass = G4std::max(minMass, G4std::min(maxMass, poleMass));
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





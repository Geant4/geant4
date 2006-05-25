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
// $Id: G4DNAProtonChargeIncrease.hh,v 1.2 2006-05-25 17:57:10 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef   G4DNAPROTONCHARGEINCREASE_HH
#define  G4DNAPROTONCHARGEINCREASE_HH 1

#include "G4DNAChargeIncreaseInWater.hh"
#include "G4DNAProtonChargeIncreaseTotalCrossSectionPolicy.hh"
#include "G4DNAProtonChargeIncreaseFinalStatesPolicy.hh"
#include "G4DNAStopAndKillBelowEnergyLimitPolicy.hh"

class G4DNAProtonChargeIncreaseEnergyLimitsPolicy
{
protected:
  G4DNAProtonChargeIncreaseEnergyLimitsPolicy();

  const G4double     lowEnergyLimit;
  const G4bool       zeroBelowLowEnergyLimit;
  const G4double     highEnergyLimit;
  const G4bool       zeroAboveHighEnergyLimit;
};

class G4DNAProtonChargeIncreaseIncomingParticlePolicy
{
protected:
  G4DNAProtonChargeIncreaseIncomingParticlePolicy();
  const G4ParticleDefinition* IncomingParticleDefinition(void) const;
};

class G4DNAProtonChargeIncrease : public G4DNAChargeIncreaseInWater<G4DNAProtonChargeIncreaseTotalCrossSectionPolicy<G4DNAProtonChargeIncreaseIncomingParticlePolicy, G4DNAProtonChargeIncreaseEnergyLimitsPolicy>, G4DNAProtonChargeIncreaseFinalStatesPolicy<G4DNAProtonChargeIncreaseEnergyLimitsPolicy> >
{
public:
  G4DNAProtonChargeIncrease(const G4String& name = "G4DNAProtonChargeIncrease") : 
    G4DNAChargeIncreaseInWater<G4DNAProtonChargeIncreaseTotalCrossSectionPolicy<G4DNAProtonChargeIncreaseIncomingParticlePolicy, G4DNAProtonChargeIncreaseEnergyLimitsPolicy>, 
			       G4DNAProtonChargeIncreaseFinalStatesPolicy<G4DNAProtonChargeIncreaseEnergyLimitsPolicy> > 
  (name) {}

  virtual ~G4DNAProtonChargeIncrease() {}
};
#endif /* G4DNAPROTONCHARGEINCREASE_HH */

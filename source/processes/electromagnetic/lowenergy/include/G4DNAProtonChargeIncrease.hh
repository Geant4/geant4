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
// $Id: G4DNAProtonChargeIncrease.hh,v 1.3 2006-06-29 19:34:59 gunter Exp $
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

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
// $Id: G4DNAHydrogenChargeIncrease.hh,v 1.1 2005-12-20 13:40:25 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef   G4DNAHYDROGENCHARGEINCREASE_HH
 #define  G4DNAHYDROGENCHARGEINCREASE_HH 1

 #include "G4DNAChargeIncreaseInWater.hh"
 #include "G4DNAHydrogenChargeIncreaseTotalCrossSectionPolicy.hh"
 #include "G4DNAHydrogenChargeIncreaseFinalStatesPolicy.hh"
 #include "G4DNAStopAndKillBelowEnergyLimitPolicy.hh"

 class G4DNAHydrogenChargeIncreaseEnergyLimitsPolicy
 {
  protected:
                      G4DNAHydrogenChargeIncreaseEnergyLimitsPolicy();

   const G4double     lowEnergyLimit;
   const G4bool       zeroBelowLowEnergyLimit;
   const G4double     highEnergyLimit;
   const G4bool       zeroAboveHighEnergyLimit;
 };

 class G4DNAHydrogenChargeIncreaseIncomingParticlePolicy
 {
  protected:
                                        G4DNAHydrogenChargeIncreaseIncomingParticlePolicy();
   const G4ParticleDefinition *         IncomingParticleDefinition(void) const;
 };

 class G4DNAHydrogenChargeIncrease : public G4DNAChargeIncreaseInWater<G4DNAHydrogenChargeIncreaseTotalCrossSectionPolicy<G4DNAHydrogenChargeIncreaseIncomingParticlePolicy, G4DNAHydrogenChargeIncreaseEnergyLimitsPolicy>, G4DNAHydrogenChargeIncreaseFinalStatesPolicy<G4DNAHydrogenChargeIncreaseEnergyLimitsPolicy> >
 {
  public:
                                         G4DNAHydrogenChargeIncrease(const G4String & name = "G4DNAHydrogenChargeIncrease") : G4DNAChargeIncreaseInWater<G4DNAHydrogenChargeIncreaseTotalCrossSectionPolicy<G4DNAHydrogenChargeIncreaseIncomingParticlePolicy, G4DNAHydrogenChargeIncreaseEnergyLimitsPolicy>, G4DNAHydrogenChargeIncreaseFinalStatesPolicy<G4DNAHydrogenChargeIncreaseEnergyLimitsPolicy> > (name) {}
   virtual                              ~G4DNAHydrogenChargeIncrease() {}
 };
#endif /* G4DNAHYDROGENCHARGEINCREASE_HH */

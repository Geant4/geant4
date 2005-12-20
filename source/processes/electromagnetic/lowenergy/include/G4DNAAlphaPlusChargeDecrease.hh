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
// $Id: G4DNAAlphaPlusChargeDecrease.hh,v 1.1 2005-12-20 13:41:32 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef   G4DNAALPHAPLUSCHARGEDECREASE_HH
 #define  G4DNAALPHAPLUSCHARGEDECREASE_HH 1

 #include "G4DNAChargeDecreaseInWater.hh"
 #include "G4DNADingfelderChargeChangeTotalCrossSectionPolicy.hh"
 #include "G4DNAAlphaPlusChargeDecreaseFinalStatesPolicy.hh"
 #include "G4DNAStopAndKillBelowEnergyLimitPolicy.hh"

 class G4DNAAlphaPlusChargeDecreaseEnergyLimitsPolicy
 {
  protected:
                      G4DNAAlphaPlusChargeDecreaseEnergyLimitsPolicy();

   const G4double     lowEnergyLimit;
   const G4bool       zeroBelowLowEnergyLimit;
   const G4double     highEnergyLimit;
   const G4bool       zeroAboveHighEnergyLimit;
 };

 class G4DNAAlphaPlusChargeDecreaseIncomingParticlePolicy
 {
  protected:
                                        G4DNAAlphaPlusChargeDecreaseIncomingParticlePolicy();
   const G4ParticleDefinition *         IncomingParticleDefinition(void) const;

   G4int                                NumberOfPartialCrossSections(void) const;

   G4double                             f0[1];
   G4double                             a0[1];
   G4double                             a1[1];
   G4double                             b0[1];
   G4double                             b1[1];
   G4double                             c0[1];
   G4double                             d0[1];
   G4double                             x0[1];
   G4double                             x1[1];
 };

 class G4DNAAlphaPlusChargeDecrease : public G4DNAChargeDecreaseInWater<G4DNADingfelderChargeChangeTotalCrossSectionPolicy<G4DNAAlphaPlusChargeDecreaseIncomingParticlePolicy, G4DNAAlphaPlusChargeDecreaseEnergyLimitsPolicy>, G4DNAAlphaPlusChargeDecreaseFinalStatesPolicy<G4DNAAlphaPlusChargeDecreaseEnergyLimitsPolicy> >
 {
  public:
                                         G4DNAAlphaPlusChargeDecrease(const G4String & name = "G4DNAAlphaPlusChargeDecrease") : G4DNAChargeDecreaseInWater<G4DNADingfelderChargeChangeTotalCrossSectionPolicy<G4DNAAlphaPlusChargeDecreaseIncomingParticlePolicy, G4DNAAlphaPlusChargeDecreaseEnergyLimitsPolicy>, G4DNAAlphaPlusChargeDecreaseFinalStatesPolicy<G4DNAAlphaPlusChargeDecreaseEnergyLimitsPolicy> > (name) {}
   virtual                              ~G4DNAAlphaPlusChargeDecrease() {}
 };
#endif /* G4DNAALPHAPLUSCHARGEDECREASE_HH */

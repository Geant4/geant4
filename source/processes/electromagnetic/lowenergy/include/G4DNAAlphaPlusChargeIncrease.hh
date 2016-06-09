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
// $Id: G4DNAAlphaPlusChargeIncrease.hh,v 1.2 2006/06/29 19:33:18 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $

#ifndef   G4DNAALPHAPLUSCHARGEINCREASE_HH
 #define  G4DNAALPHAPLUSCHARGEINCREASE_HH 1

 #include "G4DNAChargeIncreaseInWater.hh"
 #include "G4DNADingfelderChargeChangeTotalCrossSectionPolicy.hh"
 #include "G4DNAAlphaPlusChargeIncreaseFinalStatesPolicy.hh"
 #include "G4DNAStopAndKillBelowEnergyLimitPolicy.hh"

 class G4DNAAlphaPlusChargeIncreaseEnergyLimitsPolicy
 {
  protected:
                      G4DNAAlphaPlusChargeIncreaseEnergyLimitsPolicy();

   const G4double     lowEnergyLimit;
   const G4bool       zeroBelowLowEnergyLimit;
   const G4double     highEnergyLimit;
   const G4bool       zeroAboveHighEnergyLimit;
 };

 class G4DNAAlphaPlusChargeIncreaseIncomingParticlePolicy
 {
  protected:
                                        G4DNAAlphaPlusChargeIncreaseIncomingParticlePolicy();
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

 class G4DNAAlphaPlusChargeIncrease : public G4DNAChargeIncreaseInWater<G4DNADingfelderChargeChangeTotalCrossSectionPolicy<G4DNAAlphaPlusChargeIncreaseIncomingParticlePolicy, G4DNAAlphaPlusChargeIncreaseEnergyLimitsPolicy>, G4DNAAlphaPlusChargeIncreaseFinalStatesPolicy<G4DNAAlphaPlusChargeIncreaseEnergyLimitsPolicy> >
 {
  public:
                                         G4DNAAlphaPlusChargeIncrease(const G4String & name = "G4DNAAlphaPlusChargeIncrease") : G4DNAChargeIncreaseInWater<G4DNADingfelderChargeChangeTotalCrossSectionPolicy<G4DNAAlphaPlusChargeIncreaseIncomingParticlePolicy, G4DNAAlphaPlusChargeIncreaseEnergyLimitsPolicy>, G4DNAAlphaPlusChargeIncreaseFinalStatesPolicy<G4DNAAlphaPlusChargeIncreaseEnergyLimitsPolicy> > (name) {}
   virtual                              ~G4DNAAlphaPlusChargeIncrease() {}
 };
#endif /* G4DNAALPHAPLUSCHARGEINCREASE_HH */

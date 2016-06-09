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
// $Id: G4DNAAlphaPlusPlusExcitation.hh,v 1.2 2006/06/29 19:33:32 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $

#ifndef   G4DNAALPHAPLUSPLUSEXCITATION_HH
 #define  G4DNAALPHAPLUSPLUSEXCITATION_HH 1
 
 #include "G4DNAExcitationInWater.hh"
 #include "G4DNAMillerAndGreenExcitationTotalCrossSectionPolicy.hh"
 #include "G4DNAStopAndKillBelowEnergyLimitPolicy.hh"
 
 class G4DNAAlphaPlusPlusExcitationEnergyLimitsPolicy
 {
  protected:
                      G4DNAAlphaPlusPlusExcitationEnergyLimitsPolicy();
  
   const G4double     lowEnergyLimit;
   const G4bool       zeroBelowLowEnergyLimit;
   const G4double     highEnergyLimit;
   const G4bool       zeroAboveHighEnergyLimit;
 };
 
 class G4DNAAlphaPlusPlusExcitationIncomingParticlePolicy
 {
  protected:
                                        G4DNAAlphaPlusPlusExcitationIncomingParticlePolicy();
   const G4ParticleDefinition *         IncomingParticleDefinition(void) const;

   G4double                             slaterEffectiveCharge[3];
   G4double                             sCoefficient[3];
   const G4double                       kineticEnergyCorrection;
 };
 
 class G4DNAAlphaPlusPlusExcitation : public G4DNAExcitationInWater<G4DNAMillerAndGreenExcitationTotalCrossSectionPolicy<G4DNAAlphaPlusPlusExcitationIncomingParticlePolicy, G4DNAAlphaPlusPlusExcitationEnergyLimitsPolicy>, G4DNAStopAndKillBelowEnergyLimitPolicy<G4DNAAlphaPlusPlusExcitationEnergyLimitsPolicy> >
 {
  public:
                                         G4DNAAlphaPlusPlusExcitation(const G4String & name = "G4DNAAlphaPlusPlusExcitation") : G4DNAExcitationInWater<G4DNAMillerAndGreenExcitationTotalCrossSectionPolicy<G4DNAAlphaPlusPlusExcitationIncomingParticlePolicy, G4DNAAlphaPlusPlusExcitationEnergyLimitsPolicy>, G4DNAStopAndKillBelowEnergyLimitPolicy<G4DNAAlphaPlusPlusExcitationEnergyLimitsPolicy> > (name) {}
   virtual                              ~G4DNAAlphaPlusPlusExcitation() {}
 };
#endif /* G4DNAALPHAPLUSPLUSEXCITATION_HH */

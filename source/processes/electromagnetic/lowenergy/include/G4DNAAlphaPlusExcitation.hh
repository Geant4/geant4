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
// $Id: G4DNAAlphaPlusExcitation.hh,v 1.2 2006/06/29 19:33:24 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $

#ifndef   G4DNAALPHAPLUSEXCITATION_HH
 #define  G4DNAALPHAPLUSEXCITATION_HH 1
 
 #include "G4DNAExcitationInWater.hh"
 #include "G4DNAMillerAndGreenExcitationTotalCrossSectionPolicy.hh"
 #include "G4DNAStopAndKillBelowEnergyLimitPolicy.hh"
 
 class G4DNAAlphaPlusExcitationEnergyLimitsPolicy
 {
  protected:
                      G4DNAAlphaPlusExcitationEnergyLimitsPolicy();
  
   const G4double     lowEnergyLimit;
   const G4bool       zeroBelowLowEnergyLimit;
   const G4double     highEnergyLimit;
   const G4bool       zeroAboveHighEnergyLimit;
 };
 
 class G4DNAAlphaPlusExcitationIncomingParticlePolicy
 {
  protected:
                                        G4DNAAlphaPlusExcitationIncomingParticlePolicy();
   const G4ParticleDefinition *         IncomingParticleDefinition(void) const;

   G4double                             slaterEffectiveCharge[3];
   G4double                             sCoefficient[3];
   const G4double                       kineticEnergyCorrection;
 };
 
 class G4DNAAlphaPlusExcitation : public G4DNAExcitationInWater<G4DNAMillerAndGreenExcitationTotalCrossSectionPolicy<G4DNAAlphaPlusExcitationIncomingParticlePolicy, G4DNAAlphaPlusExcitationEnergyLimitsPolicy>, G4DNAStopAndKillBelowEnergyLimitPolicy<G4DNAAlphaPlusExcitationEnergyLimitsPolicy> >
 {
  public:
                                         G4DNAAlphaPlusExcitation(const G4String & name = "G4DNAAlphaPlusExcitation") : G4DNAExcitationInWater<G4DNAMillerAndGreenExcitationTotalCrossSectionPolicy<G4DNAAlphaPlusExcitationIncomingParticlePolicy, G4DNAAlphaPlusExcitationEnergyLimitsPolicy>, G4DNAStopAndKillBelowEnergyLimitPolicy<G4DNAAlphaPlusExcitationEnergyLimitsPolicy> > (name) {}
   virtual                              ~G4DNAAlphaPlusExcitation() {}
 };
#endif /* G4DNAALPHAPLUSEXCITATION_HH */

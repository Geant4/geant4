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
// $Id: G4DNAElectronExcitation.hh,v 1.2 2006-06-28 13:58:29 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef   G4DNAELECTRONEXCITATION_HH
 #define  G4DNAELECTRONEXCITATION_HH 1

 #include "G4DNAExcitationInWater.hh"
 #include "G4DNAEmfietzoglouExcitationTotalCrossSectionPolicy.hh"
 #include "G4DNAStopAndKillBelowEnergyLimitPolicy.hh"

 class G4DNAElectronExcitationEnergyLimitsPolicy
 {
  protected:
                      G4DNAElectronExcitationEnergyLimitsPolicy();

   const G4double     lowEnergyLimit;
   const G4bool       zeroBelowLowEnergyLimit;
   const G4double     highEnergyLimit;
   const G4bool       zeroAboveHighEnergyLimit;
 };

 class G4DNAElectronExcitationIncomingParticlePolicy
 {
  protected:
                                        G4DNAElectronExcitationIncomingParticlePolicy();
   const G4ParticleDefinition *         IncomingParticleDefinition(void) const;
 };

 class G4DNAElectronExcitation : public G4DNAExcitationInWater<G4DNAEmfietzoglouExcitationTotalCrossSectionPolicy<G4DNAElectronExcitationIncomingParticlePolicy, G4DNAElectronExcitationEnergyLimitsPolicy>, G4DNAStopAndKillBelowEnergyLimitPolicy<G4DNAElectronExcitationEnergyLimitsPolicy> >
 {
  public:
                                         G4DNAElectronExcitation(const G4String & name = "G4DNAElectronExcitation") : G4DNAExcitationInWater<G4DNAEmfietzoglouExcitationTotalCrossSectionPolicy<G4DNAElectronExcitationIncomingParticlePolicy, G4DNAElectronExcitationEnergyLimitsPolicy>, G4DNAStopAndKillBelowEnergyLimitPolicy<G4DNAElectronExcitationEnergyLimitsPolicy> > (name) {}
   virtual                              ~G4DNAElectronExcitation() {}
 };
#endif /* G4DNAELECTRONEXCITATION_HH */

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
// $Id: G4DNAElectronElasticBrenner.hh,v 1.2 2005-06-24 10:07:13 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef   G4DNAELECTRONELASTICBRENNER_HH
 #define  G4DNAELECTRONELASTICBRENNER_HH 1
 
 #include "G4DNAElectronElasticScatteringInWater.hh"
 #include "G4DNARutherfordTotalCrossSectionPolicy.hh"
 #include "G4DNABrennerAngularDistributionPolicy.hh"
 
 struct G4DNABrennerEnergyLimitsPolicy
 {
  static const G4double     lowEnergyLimit;
  static const G4bool       zeroBelowLowEnergyLimit;
  static const G4double     highEnergyLimit;
  static const G4bool       zeroAboveHighEnergyLimit;
 };
 
 class G4DNAElectronElasticBrenner : public G4DNAElectronElasticScatteringInWater<G4DNARutherfordTotalCrossSectionPolicy<G4DNABrennerEnergyLimitsPolicy>, G4DNABrennerAngularDistributionPolicy>
 {
  public:
                                         G4DNAElectronElasticBrenner(const G4String & name = "G4DNAElectronElasticBrenner") : G4DNAElectronElasticScatteringInWater<G4DNARutherfordTotalCrossSectionPolicy<G4DNABrennerEnergyLimitsPolicy>, G4DNABrennerAngularDistributionPolicy>(name) {}
   virtual                              ~G4DNAElectronElasticBrenner() {}
 };
#endif /* G4DNAELECTRONELASTICBRENNER_HH */

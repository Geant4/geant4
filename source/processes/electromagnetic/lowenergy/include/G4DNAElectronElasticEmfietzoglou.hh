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
// $Id: G4DNAElectronElasticEmfietzoglou.hh,v 1.2 2005-06-24 10:07:13 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef   G4DNAELECTRONELASTICEMFIETZOGLOU_HH
 #define  G4DNAELECTRONELASTICEMFIETZOGLOU_HH 1
 
 #include "G4DNAElectronElasticScatteringInWater.hh"
 #include "G4DNARutherfordTotalCrossSectionPolicy.hh"
 #include "G4DNAEmfietzoglouAngularDistributionPolicy.hh"
 
 struct G4DNAEmfietzoglouEnergyLimitsPolicy
 {
  static const G4double     lowEnergyLimit;
  static const G4bool       zeroBelowLowEnergyLimit;
  static const G4double     highEnergyLimit;
  static const G4bool       zeroAboveHighEnergyLimit;
 };
 
 class G4DNAElectronElasticEmfietzoglou : public G4DNAElectronElasticScatteringInWater<G4DNARutherfordTotalCrossSectionPolicy<G4DNAEmfietzoglouEnergyLimitsPolicy>, G4DNAEmfietzoglouAngularDistributionPolicy>
 {
  public:
                                         G4DNAElectronElasticEmfietzoglou(const G4String & name = "G4DNAElectronElasticEmfietzoglou") : G4DNAElectronElasticScatteringInWater<G4DNARutherfordTotalCrossSectionPolicy<G4DNAEmfietzoglouEnergyLimitsPolicy>, G4DNAEmfietzoglouAngularDistributionPolicy>(name) {}
   virtual                              ~G4DNAElectronElasticEmfietzoglou() {}
 };
#endif /* G4DNAELECTRONELASTICEMFIETZOGLOU_HH */

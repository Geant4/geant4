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
// $Id: G4DNAElectronElasticBrenner.hh,v 1.3 2006/06/29 19:33:56 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $

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

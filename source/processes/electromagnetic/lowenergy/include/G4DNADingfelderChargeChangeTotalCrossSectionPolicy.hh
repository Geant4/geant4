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
// $Id: G4DNADingfelderChargeChangeTotalCrossSectionPolicy.hh,v 1.2 2006/06/29 19:33:50 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $

#ifndef   G4DNADINGFELDERCHARGECHANGETOTALCROSSSECTIONPOLICY_HH
 #define  G4DNADINGFELDERCHARGECHANGETOTALCROSSSECTIONPOLICY_HH 1

 #include "G4DNACrossSectionDataSet.hh"

 // IncomingParticlePolicy must provide:
 //  - [protected] const G4ParticleDefinition * IncomingParticleDefinition(void);
 //  - [protected] G4int NumberOfPartialCrossSections(void)
 //  - [protected] const G4double f0[]
 //  - [protected] const G4double a0[]
 //  - [protected] const G4double a1[]
 //  - [protected] const G4double b0[]
 //  - [protected] const G4double b1[]
 //  - [protected] const G4double c0[]
 //  - [protected] const G4double d0[]
 //  - [protected] G4double x0[]
 //  - [protected] G4double x1[]

 // EnergyLimitsPolicy must provide:
 //  - [protected] const G4double lowEnergyLimit
 //  - [protected] const G4double zeroBelowLowEnergyLimit
 //  - [protected] const G4double highEnergyLimit
 //  - [protected] const G4double zeroAboveLowEnergyLimit

 template <typename IncomingParticlePolicy, typename EnergyLimitsPolicy>
 class G4DNADingfelderChargeChangeTotalCrossSectionPolicy : public IncomingParticlePolicy, public EnergyLimitsPolicy
 {
  protected:
                                        G4DNADingfelderChargeChangeTotalCrossSectionPolicy() {}
                                       ~G4DNADingfelderChargeChangeTotalCrossSectionPolicy() {}

   G4double                             TotalCrossSection(G4double k, G4int z);
   G4int                                RandomizePartialCrossSection(G4double k, G4int z);
   void                                 BuildTotalCrossSection(void) const {}
//   G4int                                StepFunction(G4double x) const;
//   G4double                             WaterBindingEnergy(G4int levelIndex) const;

  private:
   G4double                             PartialCrossSection(G4double k, G4int z, G4int index);

   // Hides default constructor and assignment operator as private
                                        G4DNADingfelderChargeChangeTotalCrossSectionPolicy(const G4DNADingfelderChargeChangeTotalCrossSectionPolicy & copy);
   G4DNADingfelderChargeChangeTotalCrossSectionPolicy & operator=(const G4DNADingfelderChargeChangeTotalCrossSectionPolicy & right);
 };

#include "G4DNADingfelderChargeChangeTotalCrossSectionPolicy.icc"
#endif /* G4DNADINGFELDERCHARGECHANGETOTALCROSSSECTIONPOLICY_HH */


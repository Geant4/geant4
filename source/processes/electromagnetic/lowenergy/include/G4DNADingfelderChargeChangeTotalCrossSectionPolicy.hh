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
// $Id: G4DNADingfelderChargeChangeTotalCrossSectionPolicy.hh,v 1.1 2005-12-20 13:41:32 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

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


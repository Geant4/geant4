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
// $Id: G4DNAProtonChargeDecreaseTotalCrossSectionPolicy.hh,v 1.1 2005-09-15 09:04:21 zfrancis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef   G4DNAPROTONCHARGEDECREASETOTALCROSSSECTIONPOLICY_HH
 #define  G4DNAPROTONCHARGEDECREASETOTALCROSSSECTIONPOLICY_HH 1

 #include "G4DNACrossSectionDataSet.hh"

 // IncomingParticlePolicy must provide:
 //  - [protected] const G4ParticleDefinition * IncomingParticleDefinition(void);

 // EnergyLimitsPolicy must provide:
 //  - [protected] const double lowEnergyLimit
 //  - [protected] const double zeroBelowLowEnergyLimit
 //  - [protected] const double highEnergyLimit
 //  - [protected] const double zeroAboveLowEnergyLimit

 template <typename IncomingParticlePolicy, typename EnergyLimitsPolicy>
 class G4DNAProtonChargeDecreaseTotalCrossSectionPolicy : public IncomingParticlePolicy, public EnergyLimitsPolicy
 {
  protected:
                                        G4DNAProtonChargeDecreaseTotalCrossSectionPolicy() {}
                                       ~G4DNAProtonChargeDecreaseTotalCrossSectionPolicy() {}

   G4double                             TotalCrossSection(G4double k, G4int z) const;
   G4int                                RandomizePartialCrossSection(G4double k) const;
   void                                 BuildTotalCrossSection(void) const {}
   G4int                                StepFunction(G4double x) const;
   G4double                             WaterBindingEnergy(G4int levelIndex) const;

   // Hides default constructor and assignment operator as private
                                        G4DNAProtonChargeDecreaseTotalCrossSectionPolicy(const G4DNAProtonChargeDecreaseTotalCrossSectionPolicy & copy);
   G4DNAProtonChargeDecreaseTotalCrossSectionPolicy & operator=(const G4DNAProtonChargeDecreaseTotalCrossSectionPolicy & right);
 };

#include "G4DNAProtonChargeDecreaseTotalCrossSectionPolicy.icc"
#endif /* G4DNAPROTONCHARGEDECREASETOTALCROSSSECTIONPOLICY_HH */


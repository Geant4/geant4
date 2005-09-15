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
// $Id: G4DNAProtonChargeIncreaseTotalCrossSectionPolicy.hh,v 1.1 2005-09-15 18:24:17 zfrancis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef   G4DNAPROTONCHARGEINCREASETOTALCROSSSECTIONPOLICY_HH
 #define  G4DNAPROTONCHARGEINCREASETOTALCROSSSECTIONPOLICY_HH 1

 #include "G4DNACrossSectionDataSet.hh"

 // IncomingParticlePolicy must provide:
 //  - [protected] const G4ParticleDefinition * IncomingParticleDefinition(void);

 // EnergyLimitsPolicy must provide:
 //  - [protected] const double lowEnergyLimit
 //  - [protected] const double zeroBelowLowEnergyLimit
 //  - [protected] const double highEnergyLimit
 //  - [protected] const double zeroAboveLowEnergyLimit

 template <typename IncomingParticlePolicy, typename EnergyLimitsPolicy>
 class G4DNAProtonChargeIncreaseTotalCrossSectionPolicy : public IncomingParticlePolicy, public EnergyLimitsPolicy
 {
  protected:
                                        G4DNAProtonChargeIncreaseTotalCrossSectionPolicy() {}
                                       ~G4DNAProtonChargeIncreaseTotalCrossSectionPolicy() {}

   G4double                             TotalCrossSection(G4double k, G4int z) const;
   void                                 BuildTotalCrossSection(void) const {}

   // Hides default constructor and assignment operator as private
                                        G4DNAProtonChargeIncreaseTotalCrossSectionPolicy(const G4DNAProtonChargeIncreaseTotalCrossSectionPolicy & copy);
   G4DNAProtonChargeIncreaseTotalCrossSectionPolicy & operator=(const G4DNAProtonChargeIncreaseTotalCrossSectionPolicy & right);
 };

#include "G4DNAProtonChargeIncreaseTotalCrossSectionPolicy.icc"
#endif /* G4DNAPROTONCHARGEINCREASETOTALCROSSSECTIONPOLICY_HH */


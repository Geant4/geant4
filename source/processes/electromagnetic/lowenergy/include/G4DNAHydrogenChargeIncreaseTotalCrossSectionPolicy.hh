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
// $Id: G4DNAHydrogenChargeIncreaseTotalCrossSectionPolicy.hh,v 1.2 2006/06/29 19:34:32 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $

#ifndef   G4DNAHYDROGENCHARGEINCREASETOTALCROSSSECTIONPOLICY_HH
 #define  G4DNAHYDROGENCHARGEINCREASETOTALCROSSSECTIONPOLICY_HH 1

 #include "G4DNACrossSectionDataSet.hh"

 // IncomingParticlePolicy must provide:
 //  - [protected] const G4ParticleDefinition * IncomingParticleDefinition(void);

 // EnergyLimitsPolicy must provide:
 //  - [protected] const double lowEnergyLimit
 //  - [protected] const double zeroBelowLowEnergyLimit
 //  - [protected] const double highEnergyLimit
 //  - [protected] const double zeroAboveLowEnergyLimit

 template <typename IncomingParticlePolicy, typename EnergyLimitsPolicy>
 class G4DNAHydrogenChargeIncreaseTotalCrossSectionPolicy : public IncomingParticlePolicy, public EnergyLimitsPolicy
 {
  protected:
                                        G4DNAHydrogenChargeIncreaseTotalCrossSectionPolicy() {}
                                       ~G4DNAHydrogenChargeIncreaseTotalCrossSectionPolicy() {}

   G4double                             TotalCrossSection(G4double k, G4int z) const;
   G4int                                RandomizePartialCrossSection(G4double k, G4int z) const;
   void                                 BuildTotalCrossSection(void) const {}

   // Hides default constructor and assignment operator as private
                                        G4DNAHydrogenChargeIncreaseTotalCrossSectionPolicy(const G4DNAHydrogenChargeIncreaseTotalCrossSectionPolicy & copy);
   G4DNAHydrogenChargeIncreaseTotalCrossSectionPolicy & operator=(const G4DNAHydrogenChargeIncreaseTotalCrossSectionPolicy & right);
 };

#include "G4DNAHydrogenChargeIncreaseTotalCrossSectionPolicy.icc"
#endif /* G4DNAHYDROGENCHARGEINCREASETOTALCROSSSECTIONPOLICY_HH */


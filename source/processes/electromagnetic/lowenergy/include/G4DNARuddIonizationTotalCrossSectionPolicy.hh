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
// $Id: G4DNARuddIonizationTotalCrossSectionPolicy, 2005/07/22 14:33:54 francis
// GEANT4 tag $Name: emlowen-V07-01-03

#ifndef  G4DNARuddIonizationTotalCrossSectionPolicy_HH
#define  G4DNARuddIonizationTotalCrossSectionPolicy_HH 1

 #include "G4DNACrossSectionDataSet.hh"

 // IncomingParticlePolicy must provide:
 //  - [protected] const G4ParticleDefinition * IncomingParticleDefinition(void);

 // EnergyLimitsPolicy must provide:
 //  - [protected] const double lowEnergyLimit
 //  - [protected] const double zeroBelowLowEnergyLimit
 //  - [protected] const double highEnergyLimit
 //  - [protected] const double zeroAboveLowEnergyLimit

 template <typename IncomingParticlePolicy, typename EnergyLimitsPolicy>
 class G4DNARuddIonizationTotalCrossSectionPolicy : public IncomingParticlePolicy, public EnergyLimitsPolicy
 {
  protected:
                                        G4DNARuddIonizationTotalCrossSectionPolicy() {}
                                       ~G4DNARuddIonizationTotalCrossSectionPolicy() {}

   // Hides default constructor and assignment operator as private
                                        G4DNARuddIonizationTotalCrossSectionPolicy(const G4DNARuddIonizationTotalCrossSectionPolicy & copy);
   G4DNARuddIonizationTotalCrossSectionPolicy & operator=(const G4DNARuddIonizationTotalCrossSectionPolicy & right);
 };

// #include "G4DNARuddIonizationTotalCrossSectionPolicy.icc"
#endif /* G4DNARuddIonizationTotalCrossSectionPolicy_HH */



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
// $Id: G4DNAEmfietzoglouExcitationTotalCrossSectionPolicy, 2005/07/26 10:39:54 francis
// GEANT4 tag $Name: emlowen-V07-01-02

#ifndef  G4DNAEmfietzoglouExcitationTotalCrossSectionPolicy_HH
#define  G4DNAEmfietzoglouExcitationTotalCrossSectionPolicy_HH 1

 #include "G4DNACrossSectionDataSet.hh"

 // IncomingParticlePolicy must provide:
 //  - [protected] const G4ParticleDefinition * IncomingParticleDefinition(void);

 // EnergyLimitsPolicy must provide:
 //  - [protected] const double lowEnergyLimit
 //  - [protected] const double zeroBelowLowEnergyLimit
 //  - [protected] const double highEnergyLimit
 //  - [protected] const double zeroAboveLowEnergyLimit

 template <typename IncomingParticlePolicy, typename EnergyLimitsPolicy>
 class G4DNAEmfietzoglouExcitationTotalCrossSectionPolicy : public IncomingParticlePolicy, public EnergyLimitsPolicy
 {
  protected:
                                        G4DNAEmfietzoglouExcitationTotalCrossSectionPolicy() {}
                                       ~G4DNAEmfietzoglouExcitationTotalCrossSectionPolicy() {}

   G4double                             TotalCrossSection(G4double k, G4int z) const;
   G4int                                RandomizePartialCrossSection(G4double k, G4int z) const;
   G4double                             EnergyConstant(G4int excitationLevelIndex) const;
   void                                 BuildTotalCrossSection(void) const {}

  private:
   G4double                             PartialCrossSection(G4double T, G4int z,G4int excitationLevelIndex) const;


   // Hides default constructor and assignment operator as private
                                        G4DNAEmfietzoglouExcitationTotalCrossSectionPolicy(const G4DNAEmfietzoglouExcitationTotalCrossSectionPolicy & copy);
   G4DNAEmfietzoglouExcitationTotalCrossSectionPolicy & operator=(const G4DNAEmfietzoglouExcitationTotalCrossSectionPolicy & right);
 };

 #include "G4DNAEmfietzoglouExcitationTotalCrossSectionPolicy.icc"
#endif /* G4DNAEmfietzoglouExcitationTotalCrossSectionPolicy_HH */


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
// $Id: G4DNAMillerAndGreenExcitationTotalCrossSectionPolicy.hh,v 1.1 2005-07-20 10:01:54 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef   G4DNAMILLERANDGREENEXCITATIONTOTALCROSSSECTIONPOLICY_HH
 #define  G4DNAMILLERANDGREENEXCITATIONTOTALCROSSSECTIONPOLICY_HH 1
 
 #include "G4DNACrossSectionDataSet.hh"
 
 // IncomingParticlePolicy must provide:
 //  - [protected] const G4ParticleDefinition * IncomingParticleDefinition(void);
 //  - [protected] const G4double slaterEffectiveCharge[3]
 //  - [protected] const G4double sCoefficient[3]
 //  - [protected] const G4double kineticEnergyCorrection;
 
 // EnergyLimitsPolicy must provide:
 //  - [protected] const double lowEnergyLimit
 //  - [protected] const double zeroBelowLowEnergyLimit
 //  - [protected] const double highEnergyLimit
 //  - [protected] const double zeroAboveLowEnergyLimit
 
 template <typename IncomingParticlePolicy, typename EnergyLimitsPolicy>
 class G4DNAMillerAndGreenExcitationTotalCrossSectionPolicy : public IncomingParticlePolicy, public EnergyLimitsPolicy
 {
  protected:
                                        G4DNAMillerAndGreenExcitationTotalCrossSectionPolicy() {}
                                       ~G4DNAMillerAndGreenExcitationTotalCrossSectionPolicy() {}
 
   G4double                             TotalCrossSection(G4double k, G4int z) const;
   G4int                                RandomizePartialCrossSection(G4double k, G4int z) const;
   G4double                             EnergyConstant(G4int excitationLevelIndex) const;
   void                                 BuildTotalCrossSection(void) const {}
   
  private:
   G4double                             PartialCrossSection(G4double t, G4int z, G4int excitationLevelIndex) const;
   
   G4double                             S_1s(G4double t, G4double energyTransferred, G4double slaterEffectiveCharge, G4double shellNumber) const;
   G4double                             S_2s(G4double t, G4double energyTransferred, G4double slaterEffectiveCharge, G4double shellNumber) const;
   G4double                             S_2p(G4double t, G4double energyTransferred, G4double slaterEffectiveCharge, G4double shellNumber) const;
   G4double                             R(G4double t, G4double energyTransferred, G4double slaterEffectiveCharge, G4double shellNumber) const;
   
   // Hides default constructor and assignment operator as private 
                                        G4DNAMillerAndGreenExcitationTotalCrossSectionPolicy(const G4DNAMillerAndGreenExcitationTotalCrossSectionPolicy & copy);
   G4DNAMillerAndGreenExcitationTotalCrossSectionPolicy & operator=(const G4DNAMillerAndGreenExcitationTotalCrossSectionPolicy & right);
 };
 
 #include "G4DNAMillerAndGreenExcitationTotalCrossSectionPolicy.icc"
#endif /* G4DNAMILLERANDGREENEXCITATIONTOTALCROSSSECTIONPOLICY_HH */


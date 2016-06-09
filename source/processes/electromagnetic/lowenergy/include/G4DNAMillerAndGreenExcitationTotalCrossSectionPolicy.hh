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
// $Id: G4DNAMillerAndGreenExcitationTotalCrossSectionPolicy.hh,v 1.2 2006/06/29 19:34:47 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $

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


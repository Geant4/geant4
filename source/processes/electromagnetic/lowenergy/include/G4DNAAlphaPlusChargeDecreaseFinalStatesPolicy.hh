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
// $Id: G4DNAAlphaPlusChargeDecreaseFinalStates, 2005/09/14 15:46:26 francis
// GEANT4 tag $Name: emlowen-V07-01-07

#ifndef  G4DNAALPHAPLUSCHARGEDECREASEFINALSTATESPOLICY_HH
#define  G4DNAALPHAPLUSCHARGEDECREASEFINALSTATESPOLICY_HH 1

 #include "G4DNACrossSectionDataSet.hh"

 // IncomingParticlePolicy must provide:
 //  - [protected] const G4ParticleDefinition * IncomingParticleDefinition(void);

 // EnergyLimitsPolicy must provide:
 //  - [protected] const double lowEnergyLimit
 //  - [protected] const double zeroBelowLowEnergyLimit
 //  - [protected] const double highEnergyLimit
 //  - [protected] const double zeroAboveLowEnergyLimit

 template <typename EnergyLimitsPolicy>
 class G4DNAAlphaPlusChargeDecreaseFinalStatesPolicy : public EnergyLimitsPolicy
 {
  protected:
                                        G4DNAAlphaPlusChargeDecreaseFinalStatesPolicy() {}
                                       ~G4DNAAlphaPlusChargeDecreaseFinalStatesPolicy() {}


  G4bool                                KillIncomingParticle(G4double energy) const;
  void                                  BuildFinalStatesData(void) const;
  G4int                                 NumberOfFinalStates(G4int finalStateIndex) const;
  G4ParticleDefinition*                 OutgoingParticleDefinition(G4int finalStateIndex);
  G4double                              WaterBindingEnergyConstant(G4int finalStateIndex) const;
  G4double                              OutgoingParticleBindingEnergyConstant(G4int finalStateIndex) const;
  
  // Hides default constructor and assignment operator as private
                                        G4DNAAlphaPlusChargeDecreaseFinalStatesPolicy(const G4DNAAlphaPlusChargeDecreaseFinalStatesPolicy & copy);
   G4DNAAlphaPlusChargeDecreaseFinalStatesPolicy & operator=(const G4DNAAlphaPlusChargeDecreaseFinalStatesPolicy & right);
 };

 #include "G4DNAAlphaPlusChargeDecreaseFinalStatesPolicy.icc"
#endif /* G4DNAALPHAPLUSCHARGEDECREASEFINALSTATESPOLICY_HH */



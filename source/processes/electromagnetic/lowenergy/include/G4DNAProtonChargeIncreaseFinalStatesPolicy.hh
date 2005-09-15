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
// $Id: G4DNAProtonChargeIncreaseFinalStates, 2005/09/15 14:20:26 Ziad FRANCIS
// GEANT4 tag $Name: emlowen-V07-01-08

#ifndef  G4DNAProtonChargeIncreaseFinalStatesPolicy_HH
#define  G4DNAProtonChargeIncreaseFinalStatesPolicy_HH 1

 #include "G4DNACrossSectionDataSet.hh"

 // IncomingParticlePolicy must provide:
 //  - [protected] const G4ParticleDefinition * IncomingParticleDefinition(void);

 // EnergyLimitsPolicy must provide:
 //  - [protected] const double lowEnergyLimit
 //  - [protected] const double zeroBelowLowEnergyLimit
 //  - [protected] const double highEnergyLimit
 //  - [protected] const double zeroAboveLowEnergyLimit

 template <typename EnergyLimitsPolicy>
 class G4DNAProtonChargeIncreaseFinalStatesPolicy : public EnergyLimitsPolicy
 {
  protected:
                                        G4DNAProtonChargeIncreaseFinalStatesPolicy() {}
                                       ~G4DNAProtonChargeIncreaseFinalStatesPolicy() {}


  G4bool                                KillIncomingParticle(G4double energy) const;
  void                                  BuildFinalStatesData(void) const;
  G4int                                 NumberOfFinalStates(void) const;
  G4ParticleDefinition*                 OutgoingParticleDefinition(G4int finalStateIndex);
  G4double                              ElectronsBindingEnergyConstant(G4int finalStateIndex) const;

  // Hides default constructor and assignment operator as private
                                        G4DNAProtonChargeIncreaseFinalStatesPolicy(const G4DNAProtonChargeIncreaseFinalStatesPolicy & copy);
   G4DNAProtonChargeIncreaseFinalStatesPolicy & operator=(const G4DNAProtonChargeIncreaseFinalStatesPolicy & right);
 };

 #include "G4DNAProtonChargeIncreaseFinalStatesPolicy.icc"
#endif /* G4DNAProtonChargeIncreaseFinalStatesPolicy_HH */


